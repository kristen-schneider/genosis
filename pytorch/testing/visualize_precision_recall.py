import matplotlib

matplotlib.use("Agg")
import argparse
from collections import defaultdict

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
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
        <sample>\t<precision>\t<recall>
    """
    data = {}
    print(type(files[0]))
    for file in files:
        with open(file, "r") as f:
            precisions = defaultdict(list)
            recalls = defaultdict(list)
            segment = int(file.split(".")[-2])
            for line in f:
                _, precision, recall = line.strip().split("\t")
                precisions[segment].append(float(precision))
                recalls[segment].append(float(recall))
            data[segment] = (np.mean(precisions[segment]), np.mean(recalls[segment]))
    return data


def plot(files: list[str], output: str):
    """
    data: dict keyed by segment; values are tuples of precision and recall.
    Make a bar plot of precision and recall by segment
    """
    data = load_data(files)
    data = pd.DataFrame.from_dict(
        data,
        orient="index",
        columns=["Precision", "Recall"],
    )

    sns.barplot(data=data)
    plt.xlabel("Segment")
    plt.ylabel("Jaccard Similarity")
    plt.title("Jaccard Similarity by Segment")
    plt.savefig(output)


if __name__ == "__main__":
    args = parse_args()
    plot(args.input, args.output)
