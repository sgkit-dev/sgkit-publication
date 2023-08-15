import matplotlib.pyplot as plt
import humanize
import pandas as pd
import click


@click.command
@click.argument("datafile", type=click.File("r"))
@click.argument("output", type=click.Path())
def data_scaling(datafile, output):
    """
    Plot the figure showing file size and (basic) processing time scaling
    with sample size.
    """
    df = pd.read_csv(datafile).sort_values("num_samples")
    print(df)

    fig, ax = plt.subplots(1, 1)
    plt.loglog(df["num_samples"], df["bcf_size"], ".-", label="bcf")
    plt.loglog(df["num_samples"], df["sgkit_size"], ".-", label="sgkit")
    plt.xlabel("Sample size (diploid)")
    plt.ylabel("File size (bytes)")

    # for n, size in zip(df["num_samples"], df["bcf_size"]):
    #     ax.annotate(
    #         f"{humanize.naturalsize(size, binary=True)}",
    #         textcoords="offset points",
    #         xytext=(-15, -15),
    #         xy=(n, 0.5),
    #         xycoords="data",
    #     )

    print(df.bcf_size / df.sgkit_size)
    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim())
    ax2.set_xscale("log")
    ax2.set_xticks(df["num_samples"])
    ax2.set_xticklabels([str(m) for m in df["num_sites"]])
    ax2.set_xlabel("Number of variants")

    ax.legend()
    plt.tight_layout()
    plt.savefig(output)


@click.group()
def cli():
    pass


cli.add_command(data_scaling)


if __name__ == "__main__":
    cli()
