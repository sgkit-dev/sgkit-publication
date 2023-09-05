import click
import sgkit.io.vcf
import tskit
import dask.distributed
import numpy as np


@click.command()
@click.argument("infile", type=click.File("rb"))
@click.argument("num_individuals", type=click.INT)
@click.argument("outfile", type=click.File("w"))
def subset_trees(infile, num_individuals, outfile):
    """
    Subset the input trees file down to the specified number of
    individuals
    """
    ts = tskit.load(infile)
    # Assuming diploids
    ts = ts.simplify(ts.samples()[: 2 * num_individuals])
    ts.dump(outfile)


# TODO replace this with a proper sgkit CLI
@click.command()
@click.argument("infile", type=click.Path())
@click.argument("outfile", type=click.Path(exists=False))
def vcf_to_sgkit_zarr(
    infile,
    outfile,
):
    """
    Convert the input from vcf/bcf to sgkit
    """
    client = dask.distributed.Client(n_workers=8, threads_per_worker=1)
    print(client)
    # Tried this, doesn't work:
    # with dask.diagnostics.ProgressBar():
    # Need temp_chunk_length for 1M samples.
    sgkit.io.vcf.vcf_to_zarr(infile, outfile, temp_chunk_length=100)
    print("done")


@click.command()
@click.argument("trees_file", type=click.Path(exists=True))
@click.argument("sg_dataset", type=click.Path(exists=True))
def verify_sgkit_conversion(trees_file, sg_dataset):
    ds = sgkit.load_dataset(sg_dataset)
    print(ds)
    ts = tskit.load(trees_file)
    assert ds.samples.shape[0] * 2 == ts.num_samples
    np.testing.assert_array_equal(ts.sites_position, ds.variant_position)
    ts_iter = ts.variants()
    ds_iter = iter(ds.call_genotype)
    with click.progressbar(range(ts.num_sites)) as bar:
        for j in bar:
            # This is probably horribly inefficient because it isn't
            # doing things chunk-wise
            row = next(ds_iter, None)
            var = next(ts_iter, None)
            G = row.to_numpy().reshape(ts.num_samples)
            np.testing.assert_array_equal(G, var.genotypes)


@click.group()
def cli():
    pass


cli.add_command(subset_trees)
cli.add_command(vcf_to_sgkit_zarr)
cli.add_command(verify_sgkit_conversion)


if __name__ == "__main__":
    cli()
