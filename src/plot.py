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
        ax.loglog(dfs["num_samples"].values, dfs["size"].values,
                ".-", color=colour, label=tool)
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
        total_cpu = dfs["user_time"].values + dfs["sys_time"].values
        ax.loglog(
            dfs["num_samples"].values,
            total_cpu,
            label=f"{tool}",
            # linestyle=ls,
            marker=".",
            color=colours[tool],
        )

        # Show wall-time too. Pipeline nature of the bcftools and genozip
        # commands means that it automatically threads, even if we don't
        # really want it to.
        ax.loglog(
            dfs["num_samples"].values,
            dfs["wall_time"].values,
            label=f"{tool}",
            linestyle=":",
            # marker=".",
            color=colours[tool],
        )

        row = dfs.iloc[-1]
        time = humanize.naturaldelta(row["user_time"] + row["sys_time"])
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
    #     for tool in colours.keys():
    #         dfs = df[(df.threads == threads) & (df.tool == tool)]
    #         # Can also plot the user-time here as a check - total usertime
    #         # should not be much affected by the number of threads
    #         ax.loglog(
    #             dfs["num_samples"],
    #             dfs["wall_time"],
    #             label=f"wall+{tool}+t{threads}",
    #             linestyle=ls,
    #             marker=".",
    #             color=colours[tool],
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
    df1 = pd.read_csv(size_data, index_col=None).sort_values("num_samples")
    df2 = pd.read_csv(time_data, index_col=False).sort_values("num_samples")

    # TODO set the width properly based on document
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(4, 6))

    plot_size(ax1, df1)
    plot_total_cpu(ax2, df2)

    ax2.set_xlabel("Sample size (diploid)")
    ax1.set_ylabel("File size (bytes)")
    ax2.set_ylabel("Time (seconds)")

    plt.tight_layout()
    plt.savefig(output)


def plot_thread_speedup(ax, df, threads):
    ax.set_title(f"{threads} threads")
    # colours = {"bcftools": bcf_colour, "sgkit": sgkit_colour, "savvy": sav_colour}
    colours = {"sgkit": sgkit_colour, "savvy": sav_colour}
    for tool in colours.keys():
        base_time = df[(df.threads == 1) & (df.tool == tool)].wall_time.values
        print(base_time)
        dfs = df[(df.threads == threads) & (df.tool == tool)]
        speedup = base_time / dfs.wall_time.values
        # Can also plot the user-time here as a check - total usertime
        # should not be much affected by the number of threads
        ax.semilogx(
            dfs["num_samples"].values,
            speedup,
            label=f"{tool}",
            # linestyle=ls,
            marker=".",
            color=colours[tool],
        )
    ax.legend()


@click.command()
@click.argument("time_data", type=click.File("r"))
@click.argument("output", type=click.Path())
def thread_speedup(time_data, output):
    df1 = pd.read_csv(time_data).sort_values("num_samples")
    print(df1)

    # # TODO set the width properly based on document
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(4, 8))

    # plot_thread_speedup(ax1, df1, 2)
    plot_thread_speedup(ax2, df1, 8)
    # plot_time(ax2, df2)

    ax2.set_xlabel("Sample size (diploid)")
    ax1.set_ylabel("Fold-speedup from 1 thread")
    ax2.set_ylabel("Fold-speedup from 1 thread")
    # ax2.set_ylabel("Time (seconds)")

    plt.tight_layout()
    plt.savefig(output)


@click.group()
def cli():
    pass


cli.add_command(data_scaling)
cli.add_command(thread_speedup)


if __name__ == "__main__":
    cli()
