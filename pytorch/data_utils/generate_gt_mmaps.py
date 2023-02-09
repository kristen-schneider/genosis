import argparse
import os

import numpy as np
from mmap_ninja.ragged import RaggedMmap


def np_generator(file: str):
    with open(file, "r") as f:
        for line in f:
            yield np.array(
                list(map(float, line.rstrip().split())),
                dtype=np.float32,
            )


def make_mmap(*, P1_file: str, P2_file: str, outP1: str, outP2: str):
    os.mkdir(outP1)
    os.mkdir(outP2)
    RaggedMmap.from_generator(
        out_dir=outP1,
        sample_generator=np_generator(P1_file),
        batch_size=1000,
        verbose=True,
    )
    RaggedMmap.from_generator(
        out_dir=outP2,
        sample_generator=np_generator(P2_file),
        batch_size=1000,
        verbose=True,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--inP1", type=str, required=True)
    parser.add_argument("--inP2", type=str, required=True)
    parser.add_argument("--outP1", type=str, required=True)
    parser.add_argument("--outP2", type=str, required=True)
    args = parser.parse_args()

    make_mmap(
        P1_file=args.inP1,
        P2_file=args.inP2,
        outP1=args.outP1,
        outP2=args.outP2,
    )
