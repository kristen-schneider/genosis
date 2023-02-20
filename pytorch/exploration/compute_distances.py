import argparse
import random

import numpy as np

N_SAMPLES = 6404  # will use this to infer shape of memmap


def parse_args():
    """
    Parse command line arguments.
    Provide a help string
    function args one per line
    """
    parser = argparse.ArgumentParser(
        description="Compute distances between samples.",
    )
    parser.add_argument(
        "--memmap",
        type=str,
        help="Memmap file.",
    )
    parser.add_argument(
        "--pairings",
        dest="pairings_file",
        type=str,
        help="Pairings file.",
    )
    parser.add_argument(
        "--output",
        type=str,
        help="Output file.",
    )
    parser.add_argument(
        "--samples_list",
        type=str,
        help="Samples list.",
    )
    return parser.parse_args()


def compute_shape(
    nrows: int,
    memmap: str,
) -> tuple[int, int]:
    """
    Compute shape of memmap.
    """
    bytes_per_row = nrows * np.dtype(np.int8).itemsize
    total_bytes = np.memmap(memmap, dtype=np.int8, mode="r").nbytes
    ncols = total_bytes // bytes_per_row
    return nrows, ncols


def load_pairings(
    pairings_path: str,
) -> dict[str, list[str]]:
    pairings = {}
    with open(pairings_path) as f:
        for line in f:
            A = line.strip().split()
            pairings[A[0]] = A[1:]
    return pairings


def cosine_similarity(
    x: np.ndarray,
    y: np.ndarray,
) -> float:
    a = x.astype(np.float32)
    b = y.astype(np.float32)
    return np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b))


def make_sample_mapping(samples: list[str]) -> dict[str, int]:
    """
    Make a mapping from sample name to index.
    """
    return {sample: i for i, sample in enumerate(samples)}


def main(
    memmap: str,
    pairings_file: str,
    output: str,
    samples_list: str,
) -> None:
    samples = [line.strip() for line in open(samples_list)]
    sample2idx = make_sample_mapping(samples)
    pairings = load_pairings(pairings_file)

    nrows, ncols = compute_shape(nrows=N_SAMPLES, memmap=memmap)
    data = np.memmap(memmap, dtype=np.int8, mode="r", shape=(nrows, ncols))

    with open(output, "w") as f:
        for sample, rest in pairings.items():
            for s in rest:
                i = sample2idx[sample]
                hap = random.choice(["_0", "_1"])
                j = sample2idx[s + hap]  # add haplotype info

                dist = cosine_similarity(
                    np.array(data[i], dtype=np.int8),
                    np.array(data[j], dtype=np.int8),
                )
                f.write(f"{sample}\t{s+hap}\t{dist}\n")


if __name__ == "__main__":
    args = parse_args()
    main(**vars(args))
