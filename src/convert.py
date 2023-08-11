import click
import sgkit.io.vcf
import tskit


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
    sgkit.io.vcf.vcf_to_zarr(infile, outfile)


@click.group()
def cli():
    pass


cli.add_command(subset_trees)
cli.add_command(vcf_to_sgkit_zarr)


if __name__ == "__main__":
    cli()
