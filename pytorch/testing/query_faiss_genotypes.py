import argparse
import os

import faiss
import numpy as np


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--queries",
        dest="queries_file",
        type=str,
        required=True,
        help="Path to the query file",
    )
    parser.add_argument(
        "--index",
        dest="index_file",
        type=str,
        required=True,
        help="Path to the index",
    )
    parser.add_argument(
        "--ids",
        dest="ids_file",
        type=str,
        required=True,
        help="Path to the ids",
    )
    parser.add_argument(
        "--output",
        type=str,
        required=True,
        help="Path of output file",
    )
    parser.add_argument(
        "--k",
        type=int,
        required=True,
        help="Number of neighbors to return",
    )
    return parser.parse_args()


def get_segment(queries_file: str):
    """
    Get the segment number from the index path
    File format:  $path_to/1KG.data.seg.{segment}.encoded
    """
    return int(queries_file.split(".")[-2])


def get_queries(query_path: str):
    """
    Get the queries from the query file
    File format: $sample_id\t$query (space separated ints)
    """

    samples = []
    queries = []
    with open(query_path, "r") as f:
        for line in f:
            A = line.rstrip().split()
            samples.append(A[0])
            queries.append(np.array([float(x) for x in A[1:]]))

    # sort by sample name
    samples, queries = zip(*sorted(zip(samples, queries)))
    return samples, np.array(queries, dtype=np.float32)


def get_id_mapping(ids_path: str):
    """
    Get the mapping from index_id (int) to sample_id (str)
    """
    with open(ids_path, "r") as f:
        return [line.strip() for line in f]


def main(queries_file: str, index_file: str, ids_file: str, output: str, k: int):
    # get the segment number
    os.makedirs(os.path.dirname(output), exist_ok=True)
    segment = get_segment(queries_file)
    id2sample = get_id_mapping(ids_file)
    samples, queries = get_queries(queries_file)

    index = faiss.read_index(index_file)
    D, I = index.search(queries, k)

    # write the results
    with open(output, "w") as f:
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
