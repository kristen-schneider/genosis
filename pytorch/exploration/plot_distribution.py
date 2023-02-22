import argparse

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

# f"""
# python exploration/plot_distribution.py \
#   --segment {segment} \
#   --distances {{input}} \
#   --output {{output}}
# """


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--segment",
        type=str,
        required=True,
        help="segment name",
    )
    parser.add_argument(
        "--distances",
        dest="distances_path",
        type=str,
        required=True,
        help="Path to the distances file",
    )
    parser.add_argument(
        "--output",
        type=str,
        required=True,
        help="Path to the output file",
    )
    return parser.parse_args()


def load_data(distances_path: str) -> np.ndarray:
    with open(distances_path, "r") as f:
        return np.array([float(line.rstrip().split()[2]) for line in f])


def plot_distribution(*, distances: np.ndarray, segment:str, output: str):
    sns.distplot(distances, kde=False)
    plt.xlabel("Distance (m)")
    plt.ylabel("Count")
    plt.title(f"Distribution of distances for segment {segment}")
    plt.savefig(output)


if __name__ == "__main__":
    args = parse_args()
    distances = load_data(args.distances_path)
    plot_distribution(
        distances=distances,
        segment=args.segment,
        output=args.output,
    )
