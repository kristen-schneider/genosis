import numpy as np
from scipy.spatial.distance import cosine

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

def main(encodings: str, segment: int, output: str):
    S, E = load_encodings(encodings, segment)
    print(f"Loaded {len(S)} samples")
    print(f"{E.shape=}")

    norm = np.linalg.norm(E[0])
    print(cosine(E[0], E[0]))
    E[0] /= norm
    print(E[0].T @ E[0])

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--encodings", type=str, required=True)
    parser.add_argument("--segment", type=int, required=True)
    parser.add_argument("--output", type=str, required=True)
    args = parser.parse_args()
    main(**vars(args))


