import argparse

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

"""
Plot the distribution of distances as lineplot with error bars
using the mean/std of the distances for each segment.
"""


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--distances",
        dest="distances_files",
        nargs="+",
        type=str,
        required=True,
        help="Path to the distance files",
    )
    parser.add_argument(
        "--output",
        type=str,
        required=True,
        help="Path to the output file",
    )
    return parser.parse_args()


def get_dist_stats(distance_file: str):
    with open(distance_file, "r") as f:
        distances = [float(line.rstrip().split()[2]) for line in f]
    return np.mean(distances), np.std(distances)


def plot_stats(
    means: list[float],
    stds: list[float],
    segments: list[str],
    output: str,
):
    """
    Plot the means as a lineplot and the stds as error bars.
    """

    plt.plot(
        range(len(means)),
        means,
        color="blue",
        linewidth=2,
    )
    plt.fill_between(
        x=range(len(means)),
        y1=[means[i] + stds[i] for i in range(len(means))],
        y2=[means[i] - stds[i] for i in range(len(means))],
        color="blue",
        alpha=0.2,
    )

    plt.xticks(range(len(segments)), segments)
    plt.xlabel("Segment")
    plt.ylabel("Distance (m)")
    plt.title("Distribution of distances for each segment")
    plt.savefig(output)


def main(distances_files: list[str], output: str):
    means = []
    stds = []
    segments = []
    for distance_file in distances_files:
        segments.append(distance_file.split(".")[-2])
        mean, std = get_dist_stats(distance_file)
        means.append(mean)
        stds.append(std)
    plot_stats(means, stds, segments, output)


if __name__ == "__main__":
    args = parse_args()
    main(args.distances_files, args.output)
