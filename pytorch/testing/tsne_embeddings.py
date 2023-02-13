import argparse
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.manifold import TSNE


def parse_args():
    parser = argparse.ArgumentParser(description="TSNE embeddings")
    parser.add_argument(
        "--embeddings",
        type=str,
        required=True,
        help="Path to the embeddings file",
    )
    parser.add_argument(
        '--output',
        type=str,
        required=True,
        help='Path to the output file'
    )
    return parser.parse_args()


# Load the data
def load_data(file):
    data = []
    segments = []
    with open(file, "r") as f:
        for line in f:
            A = line.rstrip().split()
            segments.append(int(A[1]))
            data.append(np.array(A[2:]).astype(np.float32))

    return segments, np.array(data)

def plot_data(embeddings, output):
    segments, data = load_data(embeddings)

    # Reduce the dimensionality of the embeddings
    tsne = TSNE(n_components=2, verbose=1, perplexity=40, n_iter=300)
    tsne_results = tsne.fit_transform(data)

    # Plot the embeddings
    plt.figure(figsize=(16, 10))
    plt.scatter(tsne_results[:, 0], tsne_results[:, 1], c=segments)
    plt.savefig(output)


if __name__ == "__main__":
    args = parse_args()
    plot_data(args.embeddings, args.output)
