# Experimental sequential multithreaded version of VCF to sgkit zarr
# conversion. Currently quite minimal, but could use the more general
# infrastructure in sgkit (i.e. for detecting fields, etc) to expand.

import concurrent.futures
import time
import dataclasses
from typing import Any

import click
import numpy as np
import cyvcf2
import zarr
import numcodecs
import tqdm
import threading
import queue

numcodecs.blosc.use_threads = False


# based on https://gist.github.com/everilae/9697228
class ThreadedGenerator:
    def __init__(self, iterator, queue_maxsize=0, daemon=False):
        self._iterator = iterator
        self._sentinel = object()
        self._queue = queue.Queue(maxsize=queue_maxsize)
        self._thread = threading.Thread(name="generator_thread", target=self._run)
        self._thread.daemon = daemon

    def _run(self):
        # TODO check whether this correctly handles errors in the decode
        # thread.
        try:
            for value in self._iterator:
                self._queue.put(value)

        finally:
            self._queue.put(self._sentinel)

    def __iter__(self):
        self._thread.start()
        for value in iter(self._queue.get, self._sentinel):
            yield value

        self._thread.join()


def flush_array(executor, np_buffer, zarr_array):
    """
    Flush the specified chunk aligned buffer to the specified zarr array.
    """
    assert zarr_array.shape[1:] == np_buffer.shape[1:]

    old_len = zarr_array.shape[0]
    buff_len = np_buffer.shape[0]
    new_len = old_len + buff_len
    new_shape = tuple([new_len] + list(np_buffer.shape[1:]))
    zarr_array.resize(new_shape)

    if len(np_buffer.shape) == 1:
        futures = flush_1d_array(executor, np_buffer, zarr_array)
    else:
        futures = flush_2d_array(executor, np_buffer, zarr_array)
    print(zarr_array)
    return futures


def flush_1d_array(executor, np_buffer, zarr_array):
    def flush_chunk():
        zarr_array[-np_buffer.shape[0] :] = np_buffer

    return [executor.submit(flush_chunk)]


def flush_2d_array(executor, np_buffer, zarr_array):
    # Flush each of the chunks in the second dimension separately
    chunk_width = zarr_array.chunks[1]
    zarr_array_width = zarr_array.shape[1]

    def flush_chunk(start, stop):
        zarr_array[-np_buffer.shape[0] :, start:stop] = np_buffer[:, start:stop]

    start = 0
    futures = []
    while start < zarr_array_width:
        stop = min(start + chunk_width, zarr_array_width)
        future = executor.submit(flush_chunk, start, stop)
        futures.append(future)
        start = stop

    return futures


@dataclasses.dataclass
class BufferedArray:
    np_buffer: Any
    zarr_array: Any


def convert_vcf(vcf, zarr_store):
    chunk_length = 2000
    chunk_width = 10000

    n = len(vcf.samples)
    root = zarr.group(store=zarr_store, overwrite=True)
    compressor = numcodecs.Blosc(
        cname="zstd", clevel=7, shuffle=numcodecs.Blosc.AUTOSHUFFLE
    )
    pos_array = root.empty(
        "variant_position",
        shape=(0),
        chunks=(chunk_length),
        dtype=np.int32,
        compressor=compressor,
    )
    gt_array = root.empty(
        "call_genotype",
        shape=(0, n, 2),
        chunks=(chunk_length, chunk_width),
        dtype=np.int8,
        compressor=compressor,
    )
    gt_phased_array = root.empty(
        "call_genotype_phased",
        shape=(0, n),
        chunks=(chunk_length, chunk_width),
        dtype=bool,
        compressor=compressor,
    )

    # TODO generalise this so we're allocating the buffer and the array
    # at the same time.
    pos_buffer = np.zeros((chunk_length), dtype=np.int32)
    gt_buffer = np.zeros((chunk_length, n, 2), dtype=np.int8)
    gt_phased_buffer = np.zeros((chunk_length, n), dtype=bool)

    buffered_arrays = [
        BufferedArray(pos_buffer, pos_array),
        BufferedArray(gt_buffer, gt_array),
        BufferedArray(gt_phased_buffer, gt_phased_array),
    ]

    def flush_buffers(futures):
        # Make sure previous futures have completed
        for future in concurrent.futures.as_completed(futures):
            exception = future.exception()
            if exception is not None:
                raise exception
        futures = []
        for ba in buffered_arrays:
            futures.extend(flush_array(executor, ba.np_buffer, ba.zarr_array))
            # This is important - we need to allocate a new buffer so that the
            # we can be writing to this one while the old one is being flushed
            # in the background.
            ba.np_buffer = np.zeros_like(ba.np_buffer)
        return futures

    j = 0
    chunks = 0
    futures = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as executor:
        # 1M file
        progressbar = tqdm.tqdm(total=7254858, desc="VCF rows")
        # 10 file
        # progressbar = tqdm.tqdm(total=116230, desc="VCF rows")
        # TODO this is the wrong approach here, we need to keep
        # access to the decode thread so that we can kill it
        # appropriately when an error occurs.
        for variant in ThreadedGenerator(vcf, queue_maxsize=200):
            progressbar.update()
            # Translate write this record into numpy buffers.
            pos_buffer[j] = variant.POS
            gt = variant.genotype.array()
            assert gt.shape[1] == 3
            gt_buffer[j] = gt[:, :-1]
            gt_phased_buffer[j] = gt[:, -1]

            j += 1
            if j == chunk_length:
                futures = flush_buffers(futures)
                j = 0
        # Flush the last chunk
        for ba in buffered_arrays:
            ba.np_buffer = ba.np_buffer[:j]
        flush_buffers(futures)


@click.command
@click.argument("in_path", type=click.Path())
@click.argument("out_path", type=click.Path())
def main(in_path, out_path):
    store = zarr.DirectoryStore(out_path)
    vcf = cyvcf2.VCF(in_path)
    convert_vcf(vcf, store)


if __name__ == "__main__":
    main()
