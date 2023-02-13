import argparse
import os

import faiss
import numpy as np


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--segment",
        type=int,
        required=True,
        help="segment id",
    )
    parser.add_argument(
        "--embeddings",
        dest="embeddings_path",
        type=str,
        required=True,
        help="path to encodings file",
    )
    parser.add_argument(
        "--outdir",
        type=str,
        required=True,
        help="path to output directory",
    )
    return parser.parse_args()


def main(
    *,
    segment: int,
    embeddings_path: str,
    outdir: str,
):
    """
    Given a segment and a path to the encodings, create a faiss index and save it to disk.
    Additionally, save the segment ids (integers).

    Encodings have format
    $sample_id\t$segment_id\t$embedding
    The embedding is a space-separated list of floats.
    """
    os.makedirs(outdir, exist_ok=True)

    # Load the encodings of provided segment
    # this is not efficient, but it's just a demo

    sample_ids = []
    embeddings = []
    with open(embeddings_path, "r") as f:
        for line in f:
            A = line.strip().split()
            if int(A[1]) != segment:
                continue
            sample_ids.append(A[0])
            embeddings.append([float(x) for x in A[2:]])

    # Sort by sample id:
    # Here we sort the zipped list of sample ids and embeddings.
    # This sorts by the first element of the tuple, which is the sample id.
    # Then we unpack that into a zip so that we have two lists again... fun huh?
    sample_ids, embeddings = zip(*sorted(zip(sample_ids, embeddings)))

    embeddings = np.array(embeddings, dtype="float32")

    # Create a faiss index with inner product distance
    # NOTE: to get cosine similarity, we must use 1 - distance
    index = faiss.IndexFlatIP(embeddings.shape[1])
    index.add(embeddings)

    # Save the index and the sample ids
    faiss.write_index(index, f"{outdir}/index.{segment}.faiss")
    with open(f"{outdir}/ids.{segment}.txt", "w") as f:
        for sample_id in sample_ids:
            f.write(f"{sample_id}\n")


if __name__ == "__main__":
    args = parse_args()
    main(**vars(args))
