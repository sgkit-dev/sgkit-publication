import pathlib

import humanize
import numpy as np
import pandas as pd
import tskit
import click
import sgkit as sg


def get_file_size(file):
    return file.stat().st_size


def get_dir_size(dir):
    return sum(get_file_size(f) for f in dir.glob("**/*") if f.is_file())


def du(file):
    if file.is_file():
        return get_file_size(file)
    return get_dir_size(file)


def to_GiB(num_bytes):
    return num_bytes / (1024**3)


@click.command()
@click.argument("source_pattern", type=str)
@click.argument("output", type=click.Path())
@click.option("-s", "--suffix", default="")
def process_files(source_pattern, output, suffix):
    data = []
    for ts_path in pathlib.Path().glob(source_pattern):
        ts = tskit.load(ts_path)
        click.echo(f"{ts_path} n={ts.num_individuals}, m={ts.num_sites}")
        vcf_path = ts_path.with_suffix(".bcf")
        sg_path = ts_path.with_suffix(".sgz")
        ds = sg.load_dataset(sg_path)
        assert ts.num_samples // 2 == ds.samples.shape[0]
        assert ts.num_sites == ds.variant_position.shape[0]
        assert np.array_equal(ds.variant_position, ts.tables.sites.position.astype(int))
        data.append(
            {
                "sequence_length": int(ts.sequence_length),
                "num_samples": ds.samples.shape[0],
                "num_sites": ts.num_sites,
                "tsk_size": to_GiB(du(ts_path)),
                "bcf_size": to_GiB(du(vcf_path)),
                "sgkit_size": to_GiB(du(sg_path)),
            }
        )

        df = pd.DataFrame(data)
        df.to_csv(output)
    print(df)


@click.group()
def cli():
    pass


cli.add_command(process_files)


if __name__ == "__main__":
    cli()
