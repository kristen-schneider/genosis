import argparse
import random
from collections import defaultdict

import numpy as np


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--distances",
        type=str,
        required=True,
        help="Path to distances file",
    )
    parser.add_argument(
        "--output",
        type=str,
        required=True,
        help="Path to output file",
    )
    parser.add_argument(
        "--num_bins",
        type=int,
        required=True,
        help="Number of bins to divide the space of possible distances into",
    )
    parser.add_argument(
        "--max_samples",
        type=int,
        required=True,
        help="Maximum number of samples to take from each bin",
    )

    return parser.parse_args()



def load_distances(distances_path: str) -> list[tuple[str, str, float]]:
    """
    File format: <sample1>\t<sample2>\t<distance>\n
    Load into a list of tuples containing sample1, sample2, and distance.
    """
    with open(distances_path) as f:
        distances = []
        for line in f:
            A = line.rstrip().split()
            distances.append((A[0], A[1], float(A[2])))
    random.shuffle(distances)
    return distances


def bins_full(bins: list[int], max_samples: int) -> bool:
    """
    Check if all bins are full.
    """
    return all(bins[i] == max_samples for i in range(len(bins)))


def bin_sampling(*, distances_path, num_bins, max_samples):
    """
    For each segment, divide the space of possible distances into bins and
    sample from each bin to get a representative sample of distances.
    """
    bin_size = 1.0 / num_bins
    bins = [0 for _ in range(num_bins)]

    sampled_distances = []
    pairwise_distances = load_distances(distances_path)
    for s1, s2, d in pairwise_distances:
        bin_idx = min(int(d // bin_size), len(bins) - 1)
        if s1 == s2:
            continue
        if bins[bin_idx] < max_samples:
            sampled_distances.append((s1, s2, d))
            bins[bin_idx] += 1
        if bins_full(bins, max_samples):
            break
    return sampled_distances


def write_distances(sampled_distances, output):
    with open(output, "w") as f:
        for s1, s2, d in sampled_distances:
            f.write(f"{s1}\t{s2}\t{d}\n")


def main():
    args = parse_args()
    sampled_distances = bin_sampling(
        distances_path=args.distances,
        num_bins=args.num_bins,
        max_samples=args.max_samples,
    )
    write_distances(sampled_distances, args.output)


if __name__ == "__main__":
    main()
