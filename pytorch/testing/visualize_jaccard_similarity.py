import matplotlib

matplotlib.use("Agg")
import argparse
from collections import defaultdict

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input",
        type=str,
        nargs="+",
        required=True,
        help="Input files",
    )
    parser.add_argument(
        "--output",
        type=str,
        required=True,
        help="Output plot",
    )
    return parser.parse_args()


def load_data(files):
    """
    data format:
        <sample>\t<similarity>
    """
    data = {}
    print(type(files[0]))
    for file in files:
        with open(file, "r") as f:
            x = []
            segment = int(file.split(".")[-2])
            for line in f:
                _, similarity = line.strip().split("\t")
                x.append(float(similarity))
            data[segment] = np.mean(x)
    return data


def plot(files: list[str], output: str):
    """
    data: dict keyed by sample name, values are lists of jaccard similarities
    make a violin plot where x axis is sample name and y axis is jaccard similarity
    """
    data = load_data(files)
    sns.barplot(x=list(data.keys()), y=list(data.values()))
    plt.xlabel("Segment")
    plt.ylabel("Jaccard Similarity")
    plt.title("Jaccard Similarity by Segment")
    plt.savefig(output)


if __name__ == "__main__":
    args = parse_args()
    plot(args.input, args.output)
