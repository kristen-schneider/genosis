import argparse
import os

import faiss
import numpy as np


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--index",
        type=str,
        required=True,
        help="Path to the faiss index",
        dest="index_path",
    )
    parser.add_argument(
        "--segment",
        type=int,
        required=True,
        help="segment number",
    )
    parser.add_argument(
        "--encodings",
        type=str,
        required=True,
        help="Encoding data file",
    )
    parser.add_argument(
        "--sample-list",
        type=str,
        required=True,
        help="Sample list file for corresponding faiss index (used to map index to sample)",
    )
    parser.add_argument(
        "--k",
        type=int,
        default=10,
        help="Number of nearest neighbors to return",
    )
    parser.add_argument(
        "--outdir",
        type=str,
        required=True,
        help="Output directory path.",
    )
    return parser.parse_args()


def load_encodings(encodings: str, segment: int):
    # TODO
    S = []  # sample names
    E = []  # encodings

    with open(encodings, "r") as f:
        for line in f:
            A = line.strip().split()
            if int(A[1]) != segment:
                continue
            S.append(A[0])
            E.append([float(x) for x in A[2:]])

    # sort by sample name
    # S, E = zip(*sorted(zip(S, E)))
    return S, np.array(E, dtype=np.float32)


# NOTE: if the l2 faiss index doesn't work in comparison to the dot product faiss index,
# just use the metric that I used when computing the ground truth for training
# which was a a jaccard similarity esque metric.
def main(
    index_path: str,
    outdir: str,
    segment: int,
    encodings: str,
    sample_list: str,
    k: int,
):
    os.makedirs(outdir, exist_ok=True)

    index = faiss.read_index(index_path)
    S, E = load_encodings(encodings, segment)

    with open(sample_list, "r") as f:
        index2sample = [x.strip() for x in f]

    # query faiss index
    print(f"{type(k)=}")
    D, I = index.search(E, k=20)

    # write results to file TODO
    with open(f"{outdir}/encoding_queries_{segment}.txt", "w") as f:
        f.write(f"Query\t" + "\t".join([f"NN{i}" for i in range(k)]) + "\n")
        # NOTE: the nearest neighbor of the query will be the query itself
        # since it is also in the index.  I'll keep it in for now, but later
        # I'll skip it when reporting results.
        print(f"{I.shape=}, {D.shape=}")
        for i in range(I.shape[0]):
            query = S[i]
            results = "\t".join(
                [f"{index2sample[I[i][j]]}:{D[i][j]}" for j in range(I.shape[1])]
            )
            f.write(f"{query}\t{results}\n")


if __name__ == "__main__":
    main(**vars(parse_args()))
