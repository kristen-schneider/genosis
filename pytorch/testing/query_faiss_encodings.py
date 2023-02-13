# rule QueryEmbeddingIndex:
#   """
#   Query the embedding index for the nearest neighbors
#   """
#   input:
#     index = f"{config.outdir}/faiss_encoded/index.{{segment}}.faiss",
#     ids = f"{config.outdir}/faiss_encoded/ids.{{segment}}.txt",
#   output:
#     f"{config.outdir}/embedding_queries/queries.{{segment}}.txt"
#   threads:
#     1
#   conda:
#     "envs/faiss.yaml"
#   shell:
#     f"""
#     python testing/query_faiss.py \\
#     --index {{input.index}} \\
#     --ids {{input.ids}} \\
#     --outdir {config.outdir}/embedding_queries \\
#     --k {config.num_neighbors}
#     """

import argparse
import os

import faiss
import numpy as np


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--encodings",
        dest="encodings_path",
        type=str,
        required=True,
        help="Encoding data file used as queries",
    )
    parser.add_argument(
        "--index",
        dest="index_path",
        type=str,
        required=True,
        help="Path to the faiss index",
    )
    parser.add_argument(
        "--ids",
        dest="ids_path",
        type=str,
        required=True,
        help="Path to the ids file",
    )
    parser.add_argument(
        "--output",
        type=str,
        required=True,
        help="Path to the output file",
    )
    parser.add_argument(
        "--k",
        type=int,
        required=True,
        help="Number of nearest neighbors to return",
    )
    return parser.parse_args()


def get_segment(index_path: str):
    """
    Get the segment number from the index path
    File format:  $path_to/1KG.data.seg.{segment}.encoded
    """
    return int(index_path.split(".")[-2])


def get_queries(encodings: str, segment: int):
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
    S, E = zip(*sorted(zip(S, E)))
    return S, np.array(E, dtype=np.float32)


def get_id_mapping(ids_path: str):
    """
    Get the mapping from index_id (int) to sample_id (str)
    """
    with open(ids_path, "r") as f:
        return [line.strip() for line in f]


def main(
    index_path: str,
    ids_path: str,
    encodings_path: str,
    output: str,
    k: int,
):
    os.makedirs(os.path.dirname(output), exist_ok=True)
    segment = get_segment(index_path)
    id2sample = get_id_mapping(ids_path)
    samples, encodings = get_queries(encodings_path, segment)

    index = faiss.read_index(index_path)
    D, I = index.search(encodings, k)

    with open(output, "w") as f:
        # NOTE: the nearest neighbor of the query will be the query itself
        # since it is also in the index.  I'll keep it in for now, but later
        # I'll skip it when reporting results.
        print(f"{I.shape=}, {D.shape=}")
        for i in range(I.shape[0]):
            query = samples[i]
            results = "\t".join(
                [f"{id2sample[I[i][j]]}:{D[i][j]}" for j in range(I.shape[1])]
            )
            f.write(f"{query}\t{results}\n")

if __name__ == "__main__":
    args = parse_args()
    main(**vars(args))
