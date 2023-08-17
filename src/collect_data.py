import pathlib
import subprocess
import time
import dataclasses
import multiprocessing

import psutil
import humanize
import numpy as np
import pandas as pd
import tskit
import click
import sgkit as sg
import dask.distributed


def get_file_size(file):
    return file.stat().st_size


def get_dir_size(dir):
    return sum(get_file_size(f) for f in dir.glob("**/*") if f.is_file())


def du(file):
    if file.is_file():
        return get_file_size(file)
    return get_dir_size(file)


@dataclasses.dataclass
class ProcessTimeResult:
    wall: float
    system: float
    user: float


def run_bcftools_afdist(path, num_threads, debug=False):
    # S: Total number of CPU-seconds used by the system on behalf of
    #    the process (in kernel mode), in seconds.
    # U: Total number of CPU-seconds that the process used directly
    #    (in user mode), in seconds.
    cmd = (
        "BCFTOOLS_PLUGINS=software/bcftools-1.18/plugins "
        "/usr/bin/time -f'%S %U' "
        f"software/bcftools +af-dist --threads {num_threads} "
    )
    before = time.time()
    out = subprocess.run(cmd + str(path), shell=True, check=True, capture_output=True)
    wall_time = time.time() - before
    time_str = out.stderr.decode()
    sys_time, user_time = map(float, time_str.split())
    if debug:
        print(out.stdout.decode())
    return ProcessTimeResult(wall_time, sys_time, user_time)


def get_prob_dist(ds, num_bins=10):
    bins = np.linspace(0, 1.0, num_bins + 1)
    bins[-1] += 0.01
    af = ds.variant_allele_frequency.values[:, 1]
    pRA = 2 * af * (1 - af)
    pAA = af * af
    hets = ds.variant_n_het.values
    homs = ds.variant_n_hom_alt.values
    a = np.bincount(np.digitize(pRA, bins), weights=hets, minlength=num_bins + 1)
    b = np.bincount(np.digitize(pAA, bins), weights=homs, minlength=num_bins + 1)

    count = (a + b).astype(int)
    return pd.DataFrame({"start": bins[:-1], "stop": bins[1:], "prob_dist": count[1:]})


def sgkit_afdist_worker(ds_path, num_threads, debug, conn):
    before = time.time()
    with dask.distributed.Client(
        processes=False, threads_per_worker=num_threads
    ) as client:
        print(client)
        ds = sg.load_dataset(ds_path)
        ds = sg.variant_stats(ds)
        df = get_prob_dist(ds)
    wall_time = time.time() - before
    cpu_times = psutil.Process().cpu_times()
    print(cpu_times)
    if debug:
        print(df)
    conn.send(f"{wall_time} {cpu_times.user} {cpu_times.system}")


def run_sgkit_afdist(ds_path, num_threads, debug=False):
    conn1, conn2 = multiprocessing.Pipe()
    p = multiprocessing.Process(
        target=sgkit_afdist_worker, args=(ds_path, num_threads, debug, conn2)
    )
    p.start()
    value = conn1.recv()
    wall_time, user_time, sys_time = map(float, value.split())
    p.join()
    if p.exitcode != 0:
        raise ValueError()
    p.close()
    return ProcessTimeResult(wall_time, sys_time, user_time)


@click.command()
@click.argument("source_pattern", type=str)
@click.argument("output", type=click.Path())
@click.option("-s", "--suffix", default="")
@click.option("--debug", is_flag=True)
def processing_time(source_pattern, output, suffix, debug):
    data = []
    for ts_path in pathlib.Path().glob(source_pattern):
        ts = tskit.load(ts_path)
        click.echo(f"{ts_path} n={ts.num_individuals}, m={ts.num_sites}")
        bcf_path = ts_path.with_suffix(".bcf")
        bcf_af_path = ts_path.with_suffix(".tags.bcf")
        sg_path = ts_path.with_suffix(".sgz")
        if not sg_path.exists:
            print("Skipping missing", sg_path)
            continue
        ds = sg.load_dataset(sg_path)
        assert ts.num_samples // 2 == ds.samples.shape[0]
        assert ts.num_sites == ds.variant_position.shape[0]
        assert np.array_equal(ds.variant_position, ts.tables.sites.position.astype(int))
        for num_threads in [1, 2, 8]:
            bcf_time = run_bcftools_afdist(bcf_af_path, num_threads, debug)
            print("BCF:", bcf_time)
            sgkit_time = run_sgkit_afdist(sg_path, num_threads, debug)
            print("SG:", sgkit_time)
            for result, prog in [(bcf_time, "bcftools"), (sgkit_time, "sgkit")]:
                data.append(
                    {
                        "num_samples": ds.samples.shape[0],
                        "num_sites": ts.num_sites,
                        "prog": prog,
                        "threads": num_threads,
                        "user_time": result.user,
                        "sys_time": result.system,
                        "wall_time": result.wall,
                    }
                )

        df = pd.DataFrame(data).sort_values("num_samples")
        df.to_csv(output)
        print(df)


@click.command()
@click.argument("source_pattern", type=str)
@click.argument("output", type=click.Path())
@click.option("-s", "--suffix", default="")
@click.option("--debug", is_flag=True)
def file_size(source_pattern, output, suffix, debug):
    data = []
    for ts_path in pathlib.Path().glob(source_pattern):
        ts = tskit.load(ts_path)
        click.echo(f"{ts_path} n={ts.num_individuals}, m={ts.num_sites}")
        bcf_path = ts_path.with_suffix(".bcf")
        bcf_af_path = ts_path.with_suffix(".tags.bcf")
        sg_path = ts_path.with_suffix(".sgz")
        if not sg_path.exists:
            print("Skipping missing", sg_path)
            continue
        ds = sg.load_dataset(sg_path)
        assert ts.num_samples // 2 == ds.samples.shape[0]
        assert ts.num_sites == ds.variant_position.shape[0]
        assert np.array_equal(ds.variant_position, ts.tables.sites.position.astype(int))
        data.append(
            {
                "sequence_length": int(ts.sequence_length),
                "num_samples": ds.samples.shape[0],
                "num_sites": ts.num_sites,
                "tsk_size": du(ts_path),
                "bcf_size": du(bcf_path),
                "sgkit_size": du(sg_path),
            }
        )
        df = pd.DataFrame(data).sort_values("num_samples")
        df.to_csv(output)
        # print(df)

    print(df)


@click.command()
def tmp():
    ds = sg.load_dataset("scaling/data/chr21_n10.sgz")
    print(ds)
    ds = sg.variant_stats(ds)
    print(get_prob_dist(ds))


@click.group()
def cli():
    pass


cli.add_command(tmp)
cli.add_command(file_size)
cli.add_command(processing_time)


if __name__ == "__main__":
    cli()
