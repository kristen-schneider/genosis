import argparse
import os

import faiss
import numpy as np


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--file",
        type=str,
        required=True,
        help="genotype segment file",
    )
    parser.add_argument(
        "--outdir",
        type=str,
        required=True,
        help="output directory for index files",
    )
    return parser.parse_args()


def main(*, file, outdir):
    os.makedirs(outdir, exist_ok=True)

    sample_ids = []
    segment = file.split('.')[-2]
    genotypes = []
    with open(file) as f:
        for line in f:
            A = line.rstrip().split()
            sample_ids.append(A[0])
            genotypes.append(np.array(A[1:], dtype=np.float32))

    # sort by sample id
    sample_ids, genotypes = zip(*sorted(zip(sample_ids, genotypes)))
    genotypes = np.array(genotypes)

    # the other index we will compare with will use IP as the distance metric
    # so we this will just be for seeing relative matching performance
    # using stuff like set membership of gold standard queries vs. encoded queries
    index = faiss.IndexFlatL2(genotypes.shape[1])
    index.add(genotypes)

    faiss.write_index(index, f"{outdir}/index.{segment}.faiss")
    with open(f"{outdir}/ids.{segment}.txt", "w") as f:
        for sample_id in sample_ids:
            f.write(f"{sample_id}\n")

if __name__ == "__main__":
    main(**vars(parse_args()))
