import argparse

import numpy as np

"""
From the resampled distances, make a dataset set of pairs
of samples, over all segments.
"""


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--subtract_segment_from_pos",
        action="store_true",
        help="Subtract segment number from position",
    )
    parser.add_argument(
        "--pos_files",
        nargs="+",
        required=True,
        help="Position files",
    )
    parser.add_argument(
        "--distance_files",
        nargs="+",
        required=True,
        help="Distance file",
    )
    parser.add_argument(
        "--P1",
        required=True,
        help="Output file for P1",
    )
    parser.add_argument(
        "--P2",
        required=True,
        help="Output file for P2",
    )
    parser.add_argument(
        "--D",
        required=True,
        help="Output file for D",
    )
    return parser.parse_args()


def load_distances(distance_file: str) -> list[tuple[str, str, str]]:
    """
    For a given segment, load the list of pairwise distances with format
    sample1 sample2 distance
    """
    distances = []
    with open(distance_file, "r") as f:
        for line in f:
            seg1, seg2, distance = line.rstrip().split()
            distances.append((seg1, seg2, distance))
    return distances


def load_positions(pos_file: str, subtract_segment_from_pos: bool) -> dict[str, list[float]]:
    """
    For a given segment, load the list of position vectors with format
    sample p1 p2 ... pn
    """
    positions = {}
    segment = int(pos_file.split(".")[-2])
    with open(pos_file, "r") as f:
        if subtract_segment_from_pos:
            for line in f:
                seg, *pos = line.rstrip().split()
                positions[seg] = [float(p) - segment for p in pos]
        else:
            for line in f:
                seg, *pos = line.rstrip().split()
                positions[seg] = [float(p) for p in pos]
    return positions


def main(
    *,
    subtract_segment_from_pos: bool,
    pos_files: list[str],
    distance_files: list[str],
    P1: str,
    P2: str,
    D: str,
):
    """
    For each segment, load the list of pairwise distances and the list of
    position vectors. Then, for each pair of samples, write the position
    vectors and the distance to the output files.
    """

    # sort files by segment number
    pos_files = sorted(
        pos_files,
        key=lambda x: int(x.split(".")[-2]),
    )
    distance_files = sorted(
        distance_files,
        key=lambda x: int(x.split(".")[-2]),
    )

    with open(P1, "w") as fP1, open(P2, "w") as fP2, open(D, "w") as fD:
        for pos_file, distance_file in zip(pos_files, distance_files):
            positions = load_positions(pos_file, subtract_segment_from_pos)
            distances = load_distances(distance_file)

            for s1, s2, d in distances:
                p1 = positions[s1]
                p2 = positions[s2]
                fP1.write(" ".join(map(str,p1)) + "\n")
                fP2.write(" ".join(map(str,p2)) + "\n")
                fD.write(d + "\n")


if __name__ == "__main__":
    args = parse_args()
    main(**vars(args))
