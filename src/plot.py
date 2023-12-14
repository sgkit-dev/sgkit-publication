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

    GB = 2**30

    for tool, colour in colour_map.items():
        dfs = df[df.tool == tool]
        dfs = dfs.sort_values("num_samples")
        ax.loglog(
            dfs["num_samples"].values,
            dfs["size"].values,
            ".-",
            color=colour,
            label=tool,
        )
        row = dfs.iloc[-1]
        size = row["size"] / GB
        ax.annotate(
            f"{size:.0f}G",
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
    ax2.set_xlabel("Number of variants")
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

        hours = total_cpu[-1] // 3600

        ax.annotate(
            f"{hours:.0f}h",
            textcoords="offset points",
            xytext=(15, 0),
            xy=(row.num_samples, total_cpu[-1]),
            xycoords="data",
        )


def plot_thread_speedup(ax, df, threads):
    colours = {
        "bcftools": bcf_colour,
        "sgkit": sgkit_colour,
        "savvy": sav_colour,
        "genozip": genozip_colour,
    }
    # colours = {"sgkit": sgkit_colour, "savvy": sav_colour}
    for tool in colours.keys():
        base_time = df[(df.threads == 1) & (df.tool == tool)].wall_time.values
        dfs = df[(df.threads == threads) & (df.tool == tool)]
        speedup = base_time[: len(dfs)] / dfs.wall_time.values
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
        row = dfs.iloc[-1]
        # print(tool, "n=", row.num_samples, "wall time:", row.wall_time)
    # ax.legend()


@click.command()
@click.argument("size_data", type=click.File("r"))
@click.argument("output", type=click.Path())
def data_scaling(size_data, output):
    """
    Plot the figure showing file size.
    """
    df1 = pd.read_csv(size_data, index_col=None).sort_values("num_samples")

    # TODO set the width properly based on document
    fig, ax1 = plt.subplots(1, 1, figsize=(4, 3))

    plot_size(ax1, df1)
    # plot_total_cpu(ax2, df2)
    # plot_thread_speedup(ax3, df2, 8)

    ax1.set_xlabel("Sample size (diploid)")
    ax1.set_ylabel("File size (bytes)")

    plt.tight_layout()
    plt.savefig(output)


@click.command()
@click.argument("time_data", type=click.File("r"))
@click.argument("output", type=click.Path())
def whole_matrix_compute(time_data, output):
    """
    Plot the figure showing compute performance on whole-matrix afdist.
    """
    df1 = pd.read_csv(time_data, index_col=False).sort_values("num_samples")

    # TODO set the width properly based on document
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(4, 6))

    plot_total_cpu(ax1, df1)
    plot_thread_speedup(ax2, df1, 8)

    ax2.set_xlabel("Sample size (diploid)")
    ax1.set_ylabel("Time (seconds)")
    ax2.set_ylabel("Fold-speedup from 1 thread")

    ax1.set_title(f"Afdist CPU time")
    ax2.set_title(f"Speedup with 8 threads")
    ax2.legend()

    plt.tight_layout()
    plt.savefig(output)


@click.group()
def cli():
    pass


cli.add_command(data_scaling)
cli.add_command(whole_matrix_compute)


if __name__ == "__main__":
    cli()
