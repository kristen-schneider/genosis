import argparse
import os
import sys
from contextlib import ExitStack
from glob import glob

import numpy as np
from tqdm import tqdm


def load_gts(file):
    """
    Load genotypes into list of np.arrays for a particular segment (skipping sample names)
    """
    x = []
    with open(file, "r") as f:
        for line in f:
            x.append(np.array(list(map(int, line.strip().split()[1:]))))  # skip sample
        return x


def load_pos(file):
    """
    Load positions into list of strings for a particular segment (skipping sample names)
    """
    x = []
    seg_num = int(file.split(".")[-2])
    with open(file, "r") as f:
        for line in f:
            # skip the sample name and just load the positions
            # we subtract the segment number to make the positions relative to the segment
            x.append(
                " ".join(
                    map(
                        str,
                        np.array(line.strip().split()[1:], dtype=np.float32) - seg_num,
                    )
                )
            )
        return x


def jaccard(g1, g2):
    """
    Compute distance between two genotype vectors using jaccard similarity.
    """
    a = np.asarray(g1).astype(np.bool_)
    b = np.asarray(g2).astype(np.bool_)
    return np.double(np.logical_and(a, b).sum()) / np.double(np.logical_or(a, b).sum())


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--indir", type=str, required=True)
    parser.add_argument("--outdir", type=str, required=True)
    parser.add_argument("--seed", type=int, default=42)
    return parser.parse_args()


def main(args):
    seed = args.seed
    np.random.seed(seed)

    indir = args.indir
    outdir = args.outdir
    os.makedirs(outdir, exist_ok=True)

    # files are formated like this 1KG.data.seg.NNN.{pos_}encoded
    gt_files = sorted(
        glob(f"{indir}/*.encoded"),
        key=lambda x: int(x.split(".")[-2]),
    )
    pos_files = sorted(
        glob(f"{indir}/*.pos_encoded"),
        key=lambda x: int(x.split(".")[-2]),
    )

    with ExitStack() as stack:
        # output files
        # G1_out = stack.enter_context(open(f"{outdir}/G1.txt", "w"))
        # G2_out = stack.enter_context(open(f"{outdir}/G2.txt", "w"))
        P1_out = stack.enter_context(open(f"{outdir}/P1.txt", "w"))
        P2_out = stack.enter_context(open(f"{outdir}/P2.txt", "w"))
        D_out = stack.enter_context(open(f"{outdir}/D.txt", "w"))

        for gt_file, pos_file in tqdm(zip(gt_files, pos_files), desc="Segments"):
            G = load_gts(gt_file)
            P = load_pos(pos_file)
            assert len(G) == len(P)

            # get two random permuations of indices
            r1 = np.random.permutation(len(G))
            r2 = np.random.permutation(len(G))

            G1 = [G[i] for i in r1]
            G2 = [G[i] for i in r2]
            P1 = [P[i] for i in r1]
            P2 = [P[i] for i in r2]

            D = [jaccard(g1, g2) for g1, g2 in zip(G1, G2)]

            # write to files
            for p1, p2, d in zip(P1, P2, D):
                P1_out.write(p1 + "\n")
                P2_out.write(p2 + "\n")
                D_out.write(str(d) + "\n")


if __name__ == "__main__":
    main(parse_args())
