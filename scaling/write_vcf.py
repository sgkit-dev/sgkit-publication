import concurrent.futures
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
class ThreadedGenerator:
    def __init__(self, iterator, queue_maxsize=0, daemon=False):
        self._iterator = iterator
        self._sentinel = object()
        self._queue = queue.Queue(maxsize=queue_maxsize)
        self._thread = threading.Thread(name="generator_thread", target=self._run)
        self._thread.daemon = daemon

    def _run(self):
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


def flush_array(executor, buffer, array):
    """
    Flush the specified chunk aligned buffer to the specified zarr array.
    """
    assert array.shape[1:] == buffer.shape[1:]
    old_len = array.shape[0]
    buff_len = buffer.shape[0]
    new_len = old_len + buff_len
    new_shape = tuple([new_len] + list(buffer.shape[1:]))
    array.resize(new_shape)

    # Flush each of the chunks separately
    chunk_width = array.chunks[1]
    array_width = array.shape[1]

    def flush_chunk(start, stop):
        array[old_len:, start:stop] = buffer[:, start:stop]

    start = 0
    futures = []
    while start < array_width:
        stop = min(start + chunk_width, array_width)
        # print("Flush slice", start, stop) print(buffer[start: stop])
        # flush_chunk(start, stop)
        future = executor.submit(flush_chunk, start, stop)
        futures.append(future)
        start = stop

    # array[old_len:] = buffer
    print(array)
    return futures


def convert_vcf(vcf, zarr_store):
    chunk_length = 2000
    chunk_width = 10000

    n = len(vcf.samples)
    gt_buffer = np.zeros((chunk_length, n, 2), dtype=np.int8)
    root = zarr.group(store=zarr_store, overwrite=True)
    compressor = numcodecs.Blosc(
        cname="zstd", clevel=7, shuffle=numcodecs.Blosc.AUTOSHUFFLE
    )
    gt_array = root.empty(
        "call_genotype",
        shape=(0, n, 2),
        chunks=(chunk_length, chunk_width),
        dtype=np.int8,
        compressor=compressor,
    )
    j = 0
    chunks = 0
    futures = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as executor:
        progressbar = tqdm.tqdm(total=7254858, desc="Read VCF")
        for variant in ThreadedGenerator(vcf, queue_maxsize=200):
            progressbar.update()
            gt = variant.genotype.array()
            # assert gt.shape[1] == 3
            gt_buffer[j] = gt[:, :-1]
            j += 1
            if j == chunk_length:
                # Make sure previous futures have completed
                for future in concurrent.futures.as_completed(futures):
                    # print(future)gg
                    pass

                before = time.time()
                futures = flush_array(executor, gt_buffer, gt_array)
                duration = time.time() - before
                print(f"Flushed in {duration:.2f} seconds")
                after = time.time()

                # Important! Do not write into the old buffer.
                gt_buffer = np.zeros_like(gt_buffer)
                j = 0
                chunks += 1
            # if chunks == 5: break

    # for row in gt_array: print(row)


@click.command
@click.argument("in_path", type=click.Path())
@click.argument("out_path", type=click.Path())
def main(in_path, out_path):
    store = zarr.DirectoryStore(out_path)
    vcf = cyvcf2.VCF(in_path)
    print("START1")
    convert_vcf(vcf, store)


if __name__ == "__main__":
    main()
