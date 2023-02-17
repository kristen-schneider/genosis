import argparse
import random
from collections import defaultdict

from pprint import pprint


def parse_args():
    args = argparse.ArgumentParser()
    args.add_argument(
        "--sample_table",
        required=True,
        type=str,
        help="1000 genomes sample table",
    )
    args.add_argument(
        "--samples_list",
        required=True,
        type=str,
        help="file containing list of samples to choose from",
    )
    args.add_argument(
        "--output",
        required=True,
        type=str,
        help="output file",
    )
    args.add_argument(
        "--N",
        required=True,
        type=int,
        help="number of samples to choose",
    )
    return args.parse_args()


def get_subpops(sample_table: str) -> dict[str, list]:
    """
    Get a dictionary of subpopulations to samples (kept as sets)
    """
    # relevant columns
    SAMPLE = 0
    SUBPOP = 3

    subpop2samples = defaultdict(list)
    with open(sample_table, "r") as f:
        f.readline()  # skip header
        for line in f:
            A = line.strip().split("\t")
            subpop2samples[A[SUBPOP]].append(A[SAMPLE])
    return subpop2samples


def get_samples(samples_list: str) -> list[str]:
    samples = []
    with open(samples_list, "r") as f:
        for line in f:
            samples.append(line.strip())
    return samples


def choose_from_subpops(
    subpop2samples: dict[str, list],
    samples: list[str],
    N: int,
    output: str,
):
    """
    Randomly choose N samples from each subpopulation for each sample
    """
    with open(output, "w") as f:
        for sample in samples:
            f.write(f"{sample}\t")
            for subpop in subpop2samples.keys():
                chosen = "\t".join(
                    random.sample(
                        subpop2samples[subpop],
                        N
                        if N < len(subpop2samples[subpop])
                        else len(subpop2samples[subpop]),
                    )
                )
                f.write(f"{chosen}\t")
            f.write("\n")


def main(
    *,
    sample_table: str,
    samples_list: str,
    output: str,
    N: int,
):
    subpop2samples = get_subpops(sample_table)
    samples = get_samples(samples_list)
    choose_from_subpops(subpop2samples, samples, N, output)


if __name__ == "__main__":
    args = parse_args()
    main(**vars(args))
