import os
import pathlib
import subprocess
import time
import dataclasses
import multiprocessing
import tempfile

# Need to do this to avoid crashes with numba :cry:
os.environ["SGKIT_DISABLE_NUMBA_CACHE"] = "1"

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


def time_cli_command(cmd, debug):
    # FIXME this doesn't look like it's capturing time
    # time spent by threads correctly

    # S: Total number of CPU-seconds used by the system on behalf of
    #    the process (in kernel mode), in seconds.
    # U: Total number of CPU-seconds that the process used directly
    #    (in user mode), in seconds.
    before = time.time()
    if debug:
        print(cmd)
    out = subprocess.run(cmd, shell=True, check=True, capture_output=True)
    wall_time = time.time() - before
    time_str = out.stderr.decode()
    sys_time, user_time = map(float, time_str.split())
    if debug:
        print(out.stdout.decode())
    return ProcessTimeResult(wall_time, sys_time, user_time)


def savvy_version():
    cmd = "./software/savvy/bin/sav --version"
    out = subprocess.run(cmd, shell=True, check=True, capture_output=True)
    return out.stdout.decode()


def bcftools_version():
    cmd = "./software/bcftools --version"
    out = subprocess.run(cmd, shell=True, check=True, capture_output=True)
    return out.stdout.decode()


def genozip_version():
    cmd = "./software/genozip --version"
    out = subprocess.run(cmd, shell=True, check=True, capture_output=True)
    return out.stdout.decode()


def sgkit_version():
    return f"sgkit version: {sg.__version__}"


def get_variant_slice_region(ds, variant_slice):
    pos = ds.variant_position[variant_slice].values
    return f"{ds.contig_id.values[0]}:{pos[0]}-{pos[-1]}"


def write_sample_names(ds, sample_slice, f):
    sample_names = ds.sample_id[sample_slice].values
    for s in sample_names:
        print(s, file=f)
    f.flush()


def run_bcftools_afdist_subset(
    path, ds, variant_slice, sample_slice, *, num_threads, debug=False
):
    region = get_variant_slice_region(ds, variant_slice)
    with tempfile.NamedTemporaryFile("w") as f:
        write_sample_names(ds, sample_slice, f.file)
        cmd = (
            "export BCFTOOLS_PLUGINS=software/bcftools-1.18/plugins; "
            "/usr/bin/time -f'%S %U' "
            # Need to run the pipeline in a subshell to make sure we're
            # timing correctly
            "sh -c '"
            # bcftools view updates the AC and AN fields by default, but
            # it doesn't update AF (which +af-dist uses). So, we use
            # -I to prevent it updating, and then call fill-tags
            f"software/bcftools view -I -r {region} -S {f.name} "
            # Output uncompressed BCF to make pipeline more efficient
            f"-Ou --threads {num_threads} {path} | "
            f"software/bcftools +fill-tags -Ou --threads {num_threads} | "
            f"software/bcftools +af-dist --threads {num_threads}"
            "'"
        )
        return time_cli_command(cmd, debug)


# NOTE! These must be called on a file that has had fill-tags run on it.
def run_bcftools_afdist(path, *, num_threads, num_sites, debug=False):
    cmd = (
        "BCFTOOLS_PLUGINS=software/bcftools-1.18/plugins "
        "/usr/bin/time -f'%S %U' "
        f"./software/bcftools +af-dist --threads {num_threads} {path}"
    )
    return time_cli_command(cmd, debug)


def run_genozip_afdist(path, *, num_threads, num_sites, debug=False):
    cmd = (
        "export BCFTOOLS_PLUGINS=software/bcftools-1.18/plugins; "
        "/usr/bin/time -f'%S %U' "
        # Need to run the pipeline in a subshell to make sure we're
        # timing correctly
        "sh -c '"
        f"software/genocat --threads {num_threads} {path} | "
        f"software/bcftools +af-dist --threads {num_threads}"
        "'"
    )
    return time_cli_command(cmd, debug)


def run_genozip_afdist_subset(
    path, ds, variant_slice, sample_slice, *, num_threads, debug=False
):
    region = get_variant_slice_region(ds, variant_slice)
    # There's no "file" option for specifying samples with genozip,
    # so have to put on the command line. Hopefully this won't
    # break things
    samples = ",".join(ds.sample_id[sample_slice].values)
    cmd = (
        "export BCFTOOLS_PLUGINS=software/bcftools-1.18/plugins; "
        "/usr/bin/time -f'%S %U' "
        # Need to run the pipeline in a subshell to make sure we're
        # timing correctly
        "sh -c '"
        # f"software/bcftools view -I -r {region} -S {f.name} "
        f"software/genocat -r {region} -s {samples} "
        f"--threads {num_threads} {path} | "
        f"software/bcftools +fill-tags -Ou --threads {num_threads} | "
        f"software/bcftools +af-dist --threads {num_threads}"
        "'"
    )
    return time_cli_command(cmd, debug)


def run_savvy_afdist(path, *, num_threads, num_sites, debug=False):
    cmd = (
        "/usr/bin/time -f'%S %U' "
        f"software/savvy-afdist/sav-afdist --threads {num_threads} "
        f"--sav-file {path} --num-variants {num_sites}"
    )
    return time_cli_command(cmd, debug)


def get_prob_dist(ds, num_bins=10):
    ds = sg.variant_stats(ds, merge=False).compute()
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


def _sgkit_afdist_subset_worker(
    ds_path, variant_slice, sample_slice, num_threads, debug, conn
):
    before = time.time()
    with dask.distributed.Client(
        processes=False, threads_per_worker=num_threads
    ) as client:
        ds = sg.load_dataset(ds_path)
        if sample_slice is not None:
            assert variant_slice is not None
            ds = ds.isel(variants=variant_slice, samples=sample_slice)
        df = get_prob_dist(ds)
    wall_time = time.time() - before
    cpu_times = psutil.Process().cpu_times()
    if debug:
        print(df)
    conn.send(f"{wall_time} {cpu_times.user} {cpu_times.system}")


def sgkit_afdist_worker(ds_path, num_threads, debug, conn):
    return _sgkit_afdist_subset_worker(ds_path, None, None, num_threads, debug, conn)


def sgkit_afdist_subset_worker(
    ds_path, variant_slice, sample_slice, num_threads, debug, conn
):
    return _sgkit_afdist_subset_worker(
        ds_path, variant_slice, sample_slice, num_threads, debug, conn
    )


def run_sgkit_afdist(ds_path, *, num_threads, num_sites, debug=False):
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


def run_sgkit_afdist_subset(
    ds_path, ds, variant_slice, sample_slice, *, num_threads, debug=False
):
    conn1, conn2 = multiprocessing.Pipe()
    p = multiprocessing.Process(
        target=sgkit_afdist_subset_worker,
        args=(ds_path, variant_slice, sample_slice, num_threads, debug, conn2),
    )
    p.start()
    value = conn1.recv()
    wall_time, user_time, sys_time = map(float, value.split())
    p.join()
    if p.exitcode != 0:
        raise ValueError()
    p.close()
    return ProcessTimeResult(wall_time, sys_time, user_time)


@dataclasses.dataclass
class Tool:
    name: str
    suffix: str
    afdist_func: None
    afdist_subset_func: None
    version_func: None


all_tools = [
    Tool("savvy", ".sav", run_savvy_afdist, None, savvy_version),
    Tool("sgkit", ".sgz", run_sgkit_afdist, run_sgkit_afdist_subset, sgkit_version),
    # Making sure we run on the output of bcftools fill-tags
    Tool(
        "bcftools",
        ".tags.bcf",
        run_bcftools_afdist,
        run_bcftools_afdist_subset,
        bcftools_version,
    ),
    Tool(
        "genozip",
        ".tags.genozip",
        run_genozip_afdist,
        run_genozip_afdist_subset,
        genozip_version,
    ),
]


@click.command()
@click.argument("src", type=click.Path(), nargs=-1)
@click.argument("output", nargs=1, type=click.Path())
@click.option("-t", "--tool", multiple=True, default=[t.name for t in all_tools])
@click.option("--num-threads", type=int, default=1)
@click.option("--debug", is_flag=True)
def processing_time(src, output, tool, num_threads, debug):
    if len(src) == 0:
        raise ValueError("Need at least one input file!")
    tool_map = {t.name: t for t in all_tools}
    tools = [tool_map[tool_name] for tool_name in tool]

    data = []
    paths = [pathlib.Path(p) for p in sorted(src)]
    for ts_path in paths:
        ts = tskit.load(ts_path)
        click.echo(f"{ts_path} n={ts.num_individuals}, m={ts.num_sites}")
        sg_path = ts_path.with_suffix(".sgz")

        if not sg_path.exists:
            print("Skipping missing", sg_path)
            continue
        ds = sg.load_dataset(sg_path)
        num_sites = ts.num_sites
        assert ts.num_samples // 2 == ds.samples.shape[0]
        assert num_sites == ds.variant_position.shape[0]
        assert np.array_equal(ds.variant_position, ts.tables.sites.position.astype(int))

        for tool in tools:
            tool_path = ts_path.with_suffix(tool.suffix)
            if debug:
                print("Running:", tool)
            result = tool.afdist_func(
                tool_path, num_threads=num_threads, num_sites=num_sites, debug=debug
            )
            data.append(
                {
                    "num_samples": ds.samples.shape[0],
                    "num_sites": ts.num_sites,
                    "tool": tool.name,
                    "threads": num_threads,
                    "user_time": result.user,
                    "sys_time": result.system,
                    "wall_time": result.wall,
                }
            )
            df = pd.DataFrame(data).sort_values(["num_samples", "tool"])
            df.to_csv(output, index=False)
            print(df)


@click.command()
@click.argument("src", type=click.Path(), nargs=-1)
@click.argument("output", type=click.Path())
@click.option("--debug", is_flag=True)
def file_size(src, output, debug):
    paths = [pathlib.Path(p) for p in sorted(src)]
    data = []
    for ts_path in paths:
        ts = tskit.load(ts_path)
        click.echo(f"{ts_path} n={ts.num_individuals}, m={ts.num_sites}")
        bcf_path = ts_path.with_suffix(".bcf")
        sg_path = ts_path.with_suffix(".sgz")
        sav_path = ts_path.with_suffix(".sav")
        genozip_path = ts_path.with_suffix(".genozip")
        if not sg_path.exists:
            print("Skipping missing", sg_path)
            continue
        ds = sg.load_dataset(sg_path)
        assert ts.num_samples // 2 == ds.samples.shape[0]
        assert ts.num_sites == ds.variant_position.shape[0]
        assert np.array_equal(ds.variant_position, ts.tables.sites.position.astype(int))
        tmap = {
            "tsk": ts_path,
            "bcf": bcf_path,
            "sgkit": sg_path,
            "sav": sav_path,
            "genozip": genozip_path,
        }
        for tool, path in tmap.items():
            size = du(path)
            data.append(
                {
                    "sequence_length": int(ts.sequence_length),
                    "num_samples": ds.samples.shape[0],
                    "num_sites": ts.num_sites,
                    "tool": tool,
                    "size": size,
                }
            )
            df = pd.DataFrame(data).sort_values(["num_samples", "tool"])
            df.to_csv(output, index=False)

    print(df)


def midslice(n, k):
    """
    Return a slice of size k from the middle of an array of size n.
    """
    n2 = n // 2
    k2 = k // 2
    return slice(n2 - k2, n2 + k2)


@click.command()
@click.argument("src", type=click.Path(), nargs=-1)
@click.argument("output", nargs=1, type=click.Path())
@click.option("-t", "--tool", multiple=True, default=[t.name for t in all_tools])
@click.option("-s", "--slice-id", multiple=True, default=["n10", "n/2"])
@click.option("--num-threads", type=int, default=1)
@click.option("--debug", is_flag=True)
def subset_processing_time(src, output, tool, slice_id, num_threads, debug):
    if len(src) == 0:
        raise ValueError("Need at least one input file!")
    tool_map = {t.name: t for t in all_tools}
    tools = [tool_map[tool_name] for tool_name in tool]

    data = []
    paths = [pathlib.Path(p) for p in sorted(src)]
    for ts_path in paths:
        ts = tskit.load(ts_path)
        click.echo(f"{ts_path} n={ts.num_individuals}, m={ts.num_sites}")
        sg_path = ts_path.with_suffix(".sgz")

        if not sg_path.exists:
            print("Skipping missing", sg_path)
            continue
        ds = sg.load_dataset(sg_path)
        num_sites = ts.num_sites
        num_samples = ds.samples.shape[0]
        assert ts.num_samples // 2 == num_samples
        assert num_sites == ds.variant_position.shape[0]
        assert np.array_equal(ds.variant_position, ts.tables.sites.position.astype(int))

        slices = {
            "n10": (midslice(num_sites, 10000), midslice(num_samples, 10)),
            "n/2": (
                midslice(num_sites, 10000),
                midslice(num_samples, num_samples // 2),
            ),
        }

        for sid in slice_id:
            variant_slice, sample_slice = slices[sid]
            for tool in tools:
                tool_path = ts_path.with_suffix(tool.suffix)
                if debug:
                    print("Running:", tool)
                result = tool.afdist_subset_func(
                    tool_path,
                    ds,
                    variant_slice,
                    sample_slice,
                    num_threads=num_threads,
                    debug=debug,
                )
                data.append(
                    {
                        "num_samples": ds.samples.shape[0],
                        "num_sites": ts.num_sites,
                        "slice": sid,
                        "tool": tool.name,
                        "threads": num_threads,
                        "user_time": result.user,
                        "sys_time": result.system,
                        "wall_time": result.wall,
                    }
                )
                df = pd.DataFrame(data).sort_values(["num_samples", "tool"])
                df.to_csv(output, index=False)
                print(df)


@click.command()
def report_versions():
    for tool in all_tools:
        print(tool.name)
        print(tool.version_func())


@click.group()
def cli():
    pass


cli.add_command(file_size)
cli.add_command(processing_time)
cli.add_command(subset_processing_time)
cli.add_command(report_versions)


if __name__ == "__main__":
    cli()
