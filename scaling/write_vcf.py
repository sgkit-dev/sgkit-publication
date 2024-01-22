# Experimental sequential multithreaded version of VCF to sgkit zarr
# conversion. Currently quite minimal, but could use the more general
# infrastructure in sgkit (i.e. for detecting fields, etc) to expand.
import concurrent.futures
import dataclasses
import subprocess
import time
import functools
import multiprocessing
import os
from typing import Any
import time

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
# Needs refactoring to allow for graceful handing of
# errors in the main thread, and in the decode thread.
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


def flush_array(executor, np_buffer, zarr_array, offset):
    """
    Flush the specified chunk aligned buffer to the specified zarr array.
    """
    assert zarr_array.shape[1:] == np_buffer.shape[1:]

    if len(np_buffer.shape) == 1:
        futures = flush_1d_array(executor, np_buffer, zarr_array, offset)
    else:
        futures = flush_2d_array(executor, np_buffer, zarr_array, offset)
    return futures


def flush_1d_array(executor, np_buffer, zarr_array, offset):
    def flush_chunk():
        zarr_array[offset : offset + np_buffer.shape[0]] = np_buffer

    return [executor.submit(flush_chunk)]


def flush_2d_array(executor, np_buffer, zarr_array, offset):
    # Flush each of the chunks in the second dimension separately

    s = slice(offset, offset + np_buffer.shape[0])

    def flush_chunk(start, stop):
        zarr_array[s, start:stop] = np_buffer[:, start:stop]

    chunk_width = zarr_array.chunks[1]
    zarr_array_width = zarr_array.shape[1]
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


def write_partition(vcf_fields, zarr_path, partition):
    # print(f"process {os.getpid()} starting")
    vcf = cyvcf2.VCF(partition.path)
    # Requires cyvcf2>=0.30.27
    num_variants = vcf.num_records

    store = zarr.DirectoryStore(zarr_path)
    root = zarr.group(store=store)
    pos_array = root["variant_position"]
    gt_array = root["call_genotype"]
    gt_phased_array = root["call_genotype_phased"]

    chunk_length = gt_array.chunks[0]
    chunk_width = gt_array.chunks[1]
    n = gt_array.shape[1]

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

    def flush_buffers(offset, futures, start=0, stop=chunk_length):
        # Make sure previous futures have completed
        for future in concurrent.futures.as_completed(futures):
            exception = future.exception()
            if exception is not None:
                raise exception

        futures = []
        for ba in buffered_arrays:
            futures.extend(
                flush_array(executor, ba.np_buffer[start:stop], ba.zarr_array, offset)
            )
            # This is important - we need to allocate a new buffer so that the
            # we can be writing to this one while the old one is being flushed
            # in the background.
            ba.np_buffer = np.zeros_like(ba.np_buffer)

        return futures

    # progressbar = tqdm.tqdm(total=num_variants, desc="variants")
    # Flushing out the chunks takes less time than reading in here in the
    # main thread, so no real point in using lots of threads.
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as executor:
        offset = partition.offset
        j = offset % chunk_length
        chunk_start = j
        futures = []

        # TODO this is the wrong approach here, we need to keep
        # access to the decode thread so that we can kill it
        # appropriately when an error occurs.
        for variant in ThreadedGenerator(vcf, queue_maxsize=200):
            # Translate this record into numpy buffers. There is some compute
            # done here, but it's not releasing the GIL, so may not be worth
            # moving to threads.
            pos_buffer[j] = variant.POS
            gt = variant.genotype.array()
            assert gt.shape[1] == 3
            gt_buffer[j] = gt[:, :-1]
            gt_phased_buffer[j] = gt[:, -1]

            j += 1
            if j == chunk_length:
                futures = flush_buffers(offset, futures, start=chunk_start)
                j = 0

                offset += chunk_length - chunk_start
                chunk_start = 0
                assert offset % chunk_length == 0

            with progress_counter.get_lock():
                # TODO reduce IPC here by incrementing less often?
                # Might not be worth the hassle
                progress_counter.value += 1

        # Flush the last chunk
        flush_buffers(offset, futures, stop=j)
    # print(f"process {os.getpid()} finishing")


@dataclasses.dataclass
class VcfChunk:
    path: str
    num_records: int
    first_position: int
    offset: int = 0


@dataclasses.dataclass
class VcfFields:
    samples: list
    # TODO other stuff like sgkit does


def scan_vcfs(paths):
    chunks = []
    vcf_fields = None
    for path in tqdm.tqdm(paths, desc="Scan"):
        vcf = cyvcf2.VCF(path)
        fields = VcfFields(samples=vcf.samples)
        if vcf_fields is None:
            vcf_fields = fields
        else:
            if fields != vcf_fields:
                raise ValueError("Incompatible VCF chunks")
        record = next(vcf)
        chunks.append(
            VcfChunk(
                path=path, num_records=count_variants(path), first_position=record.POS
            )
        )

    # Assuming these are all on the same contig for now.
    chunks.sort(key=lambda x: x.first_position)
    offset = 0
    for chunk in chunks:
        chunk.offset = offset
        offset += chunk.num_records
    return vcf_fields, chunks


def create_zarr(path, vcf_fields, partitions):
    chunk_width = 10000
    chunk_length = 2000

    n = len(vcf_fields.samples)
    m = sum(partition.num_records for partition in partitions)

    store = zarr.DirectoryStore(path)
    root = zarr.group(store=store, overwrite=True)
    compressor = numcodecs.Blosc(
        cname="zstd", clevel=7, shuffle=numcodecs.Blosc.AUTOSHUFFLE
    )
    pos_array = root.empty(
        "variant_position",
        shape=(m),
        chunks=(chunk_length),
        dtype=np.int32,
        compressor=compressor,
    )
    gt_array = root.empty(
        "call_genotype",
        shape=(m, n, 2),
        chunks=(chunk_length, chunk_width),
        dtype=np.int8,
        compressor=compressor,
    )
    gt_phased_array = root.empty(
        "call_genotype_phased",
        shape=(m, n),
        chunks=(chunk_length, chunk_width),
        dtype=bool,
        compressor=compressor,
    )


def update_bar(progress_counter, num_variants):
    pbar = tqdm.tqdm(total=num_variants)

    while (total := progress_counter.value) < num_variants:
        inc = total - pbar.n
        pbar.update(inc)
        time.sleep(0.1)


def init_workers(counter):
    global progress_counter
    progress_counter = counter


@click.command
@click.argument("vcfs", nargs=-1, required=True)
@click.argument("out_path", type=click.Path())
def main(vcfs, out_path):
    # TODO add a try-except here for KeyboardInterrupt which will kill
    # various things and clean-up.
    vcf_fields, partitions = scan_vcfs(vcfs)
    # TODO write the Zarr to a temporary name, only renaming at the end
    # on success.
    create_zarr(out_path, vcf_fields, partitions)

    total_variants = sum(partition.num_records for partition in partitions)
    global progress_counter
    progress_counter = multiprocessing.Value("i", 0)

    # start update progress bar process
    # daemon= parameter is set to True so this process won't block us upon exit
    # TODO move to thread, no need for proc
    bar_process = multiprocessing.Process(
        target=update_bar, args=(progress_counter, total_variants), name="progress"
    )  # , daemon=True)
    bar_process.start()

    f = functools.partial(write_partition, vcf_fields, out_path)

    # for partition in partitions:
    #     f(partition)

    with concurrent.futures.ProcessPoolExecutor(
        max_workers=8, initializer=init_workers, initargs=(progress_counter,)
    ) as executor:
        # Ensure we don't write collide on writing the first and last
        # chunks of partitions by doing odd and even partitions separately.

        # NOTE: this isn't a great way to do it because we tend to get gluts of
        # flushes happening at the same time. We should submit the
        # first half with slight delays, so that procs are working at offsets.
        # Waiting for all the evens to complete before starting the first odd
        # does also lead to a noticable drop in throughput halfway through.
        for parts in [partitions[::2], partitions[1::2]]:
            futures = [executor.submit(f, part) for part in parts]
            for future in concurrent.futures.as_completed(futures):
                exception = future.exception()
                if exception is not None:
                    raise exception

    assert process_counter.value == total_variants


if __name__ == "__main__":
    main()
