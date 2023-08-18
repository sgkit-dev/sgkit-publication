import matplotlib
import matplotlib.pyplot as plt
import humanize
import pandas as pd
import click

sgkit_colour = "tab:blue"
bcf_colour = "tab:orange"


def plot_size(ax, df):
    ax.loglog(df["num_samples"], df["bcf_size"], ".-", color=bcf_colour, label="bcf")
    ax.loglog(
        df["num_samples"], df["sgkit_size"], ".-", color=sgkit_colour, label="sgkit"
    )
    ax.legend()
    # for n, size in zip(df["num_samples"], df["bcf_size"]):
    #     ax.annotate(
    #         f"{humanize.naturalsize(size, binary=True)}",
    #         textcoords="offset points",
    #         xytext=(-15, -15),
    #         xy=(n, 0.5),
    #         xycoords="data",
    #     )

    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim())
    ax2.set_xscale("log")
    ax2.set_xticks(df["num_samples"])
    ax2.set_xticklabels([str(m) for m in df["num_sites"]])
    ax2.set_xlabel("Number of variants")


def plot_time(ax, df):
    colours = {"bcftools": bcf_colour, "sgkit": sgkit_colour}
    threads_lines = []
    threads_labels = []
    for threads, ls in zip([1, 2, 8], ["solid", "dashed", "dotted"]):
        threads_lines.append(
            matplotlib.lines.Line2D([0], [0], lw=2, color="black", linestyle=ls)
        )
        threads_labels.append(f"{threads} threads")
        for prog in colours.keys():
            dfs = df[(df.threads == threads) & (df.prog == prog)]
            # Can also plot the user-time here as a check - total usertime
            # should not be much affected by the number of threads
            ax.loglog(
                dfs["num_samples"],
                dfs["wall_time"],
                label=f"wall+{prog}+t{threads}",
                linestyle=ls,
                marker=".",
                color=colours[prog],
            )

    lines = [
        matplotlib.lines.Line2D([0], [0], color=colour, lw=2)
        for colour in colours.values()
    ]
    l1 = ax.legend(lines, list(colours.keys()))
    ax.add_artist(l1)
    ax.legend(threads_lines, threads_labels, loc="lower right")


@click.command()
@click.argument("size_data", type=click.File("r"))
@click.argument("time_data", type=click.File("r"))
@click.argument("output", type=click.Path())
def data_scaling(size_data, time_data, output):
    """
    Plot the figure showing file size and (basic) processing time scaling
    with sample size.
    """
    df1 = pd.read_csv(size_data).sort_values("num_samples")
    df2 = pd.read_csv(time_data).sort_values("num_samples")

    # TODO set the width properly based on document
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(4, 6))

    plot_size(ax1, df1)
    plot_time(ax2, df2)

    ax2.set_xlabel("Sample size (diploid)")
    ax1.set_ylabel("File size (bytes)")
    ax2.set_ylabel("Time (seconds)")

    plt.tight_layout()
    plt.savefig(output)


@click.group()
def cli():
    pass


cli.add_command(data_scaling)


if __name__ == "__main__":
    cli()
