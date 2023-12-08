import matplotlib
import matplotlib.pyplot as plt
import humanize
import pandas as pd
import click

sgkit_colour = "tab:blue"
bcf_colour = "tab:orange"
sav_colour = "tab:red"
genozip_colour = "tab:purple"


def plot_size(ax, df):
    colour_map = {
        "bcf": bcf_colour,
        "sgkit": sgkit_colour,
        "sav": sav_colour,
        "genozip": genozip_colour,
    }

    for tool, colour in colour_map.items():
        dfs = df[df.tool == tool]
        dfs = dfs.sort_values("num_samples")
        ax.loglog(dfs["num_samples"], dfs["size"], ".-", color=colour, label=tool)
        row = dfs.iloc[-1]
        size = humanize.naturalsize(row["size"])
        ax.annotate(
            size,
            textcoords="offset points",
            xytext=(15, 0),
            xy=(row.num_samples, row["size"]),
            xycoords="data",
        )

    num_samples = dfs["num_samples"].values
    num_sites = dfs["num_sites"].values

    ax.legend()
    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim())
    ax2.set_xscale("log")
    ax2.set_xticks(num_samples)
    ax2.set_xticklabels([humanize.metric(m) for m in num_sites])


def plot_total_cpu(ax, df):
    colours = {
        "bcftools": bcf_colour,
        "sgkit": sgkit_colour,
        "savvy": sav_colour,
        "genozip": genozip_colour,
    }

    for tool in colours.keys():
        dfs = df[(df.threads == 1) & (df.tool == tool)]
        ax.loglog(
            dfs["num_samples"],
            dfs["user_time"],
            label=f"{tool}",
            # linestyle=ls,
            marker=".",
            color=colours[tool],
        )

        row = dfs.iloc[-1]
        time = humanize.naturaldelta(row["user_time"])
        ax.annotate(
            time,
            textcoords="offset points",
            xytext=(15, 0),
            xy=(row.num_samples, row["user_time"]),
            xycoords="data",
        )

    # threads_lines = []
    # threads_labels = []
    # # for threads, ls in zip([1, 2, 8], ["solid", "dashed", "dotted"]):
    #     threads_lines.append(
    #         matplotlib.lines.Line2D([0], [0], lw=2, color="black", linestyle=ls)
    #     )
    #     threads_labels.append(f"{threads} threads")
    #     for prog in colours.keys():
    #         dfs = df[(df.threads == threads) & (df.prog == prog)]
    #         # Can also plot the user-time here as a check - total usertime
    #         # should not be much affected by the number of threads
    #         ax.loglog(
    #             dfs["num_samples"],
    #             dfs["wall_time"],
    #             label=f"wall+{prog}+t{threads}",
    #             linestyle=ls,
    #             marker=".",
    #             color=colours[prog],
    #         )

    # lines = [
    #     matplotlib.lines.Line2D([0], [0], color=colour, lw=2)
    #     for colour in colours.values()
    # ]
    # l1 = ax.legend(lines, list(colours.keys()))
    # ax.add_artist(l1)
    # ax.legend(threads_lines, threads_labels, loc="lower right")


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
    plot_total_cpu(ax2, df2)

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
