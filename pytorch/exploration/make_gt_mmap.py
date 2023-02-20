import argparse

import numpy as np


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--gts",
        type=str,
        required=True,
        help="path to genotypes",
    )
    parser.add_argument(
        "--output",
        type=str,
        required=True,
        help="path to output mmap",
    )
    return parser.parse_args()


def load_gts(gts_path: str):
    gts = []
    with open(gts_path, "r") as f:
        for line in f:
            A = line.rstrip().split()
            gts.append(np.array([int(x) for x in A[1:]], dtype=np.int8))  # skip sample
    return np.array(gts, dtype=np.int8)


def make_memmaps(arr: np.ndarray, output: str):
    """Write the gts array to a memory-mapped file."""
    gts_mmap = np.memmap(output, dtype=np.int8, mode="w+", shape=arr.shape)
    gts_mmap[:] = arr[:]
    gts_mmap.flush()


def main(gts: str, output: str):
    arr = load_gts(gts)
    make_memmaps(arr, output)


if __name__ == "__main__":
    args = parse_args()
    main(**vars(args))
