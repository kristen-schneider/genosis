import argparse
from collections import defaultdict
import matplotlib.pyplot as plt

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--faiss_enc')
    parser.add_argument('--faiss_emb')
    parser.add_argument('--plink')
    parser.add_argument('--ed')
    parser.add_argument('--train')
    parser.add_argument('--test')

    return parser.parse_args()

def main():
    args = get_args()
    print('reading samples...')
    databaseIDs = read_samples(args.train)
    queryIDs = read_samples(args.test)
    
    print('reading plink...')
    plink_sorted_scores = read_plink(args.plink, databaseIDs, queryIDs)
    print(len(plink_sorted_scores))

    print('reading euclidean distances...')
    ed_sorted_scores = read_euclidean_distance(args.ed, databaseIDs, queryIDs)
    print(len(ed_sorted_scores))
    
    print('reading faiss encoding...')
    faiss_sorted_encoding = read_faiss_encoding(args.faiss_enc)
    print(len(faiss_sorted_encoding))
    
    print('reading faiss embedding...')
    faiss_sorted_embedding = read_faiss_embedding(args.faiss_emb)
    print(len(faiss_sorted_embedding))
    x=1


def read_faiss_encoding(faiss_encoding_file):
    F = {}
    header = 9
    with open(faiss_encoding_file) as f:
        for i in range(header):
           f.readline()
        for l in f:
            if 'QUERY' in l:
                Q = l.strip().split()[1]
                continue
            try:
                A = l.rstrip().split()
                s = A[0]
                d = float(A[2])
                try:
                    F[Q].append((s, d))
                except KeyError:
                    F[Q] = [(s, d)]
            except IndexError:
                continue

        f.close()
        # sort by distance
        for faiss in F:
            F[faiss].sort(key=lambda tup: tup[1], reverse=False)
        return F

def read_faiss_embedding(faiss_embedding_file):
    F = {}
    header = 12
    with open(faiss_embedding_file) as f:
        for i in range(header):
            f.readline()
        for l in f:
            if 'QUERY' in l:
                Q = l.strip().split()[1]
                continue
            try:
                A = l.rstrip().split()
                s = A[0]
                d = float(A[2])
                try:
                    F[Q].append((s, d))
                except KeyError:
                    F[Q] = [(s, d)]
            except IndexError:
                continue

        f.close()
        # sort by distance
        for faiss in F:
            F[faiss].sort(key=lambda tup: tup[1], reverse=False)
        return F

def read_plink(plink_file, databaseIDs, queryIDs):
    header_len = 1
    P = {}
    with open(plink_file, 'r') as f:
        for l in f:
            if "time" in l:
                break;
            if header_len > 0:
                header_len = header_len - 1
                continue
            A = l.rstrip().split()
            a = A[1]
            b = A[3]
            dist = float(A[11])

            if a+'_0' in queryIDs or a+'_1' in queryIDs:
                if a not in P:
                    P[a] = []
                if b+'_0' in databaseIDs or b+'_1' in databaseIDs:
                    P[a].append((b, dist))

            if b+'_0' in queryIDs or b+'_1' in queryIDs:
                if b not in P:
                    P[b] = []
                if a+'_0' in databaseIDs or a+'_1' in databaseIDs:
                    P[b].append((a, dist))

        f.close()
        # sort by distance
        for p in P:
            P[p].sort(key=lambda tup: tup[1], reverse=True)

        return P

def read_euclidean_distance(euclidean_file, databaseIDs, queryIDs):
    E = {}
    with open(euclidean_file, 'r') as f:
        for l in f:
            if 'Query' in l:
                Q = l.strip().split()[1]
                header = f.readline()
                continue
            A = l.rstrip().split()
            s = A[0]
            d = float(A[1])
            if Q in queryIDs:
                if s in databaseIDs:
                    try:
                        E[Q].append((s, d))
                    except KeyError:
                        E[Q] = [(s, d)]
        f.close()
        # sort by distance
        for e in E:
            E[e].sort(key=lambda tup: tup[1], reverse=False)

        return E

def read_samples(samples_file):
    samples = []
    f = open(samples_file, 'r')
    for line in f:
        s = line.strip()
        samples.append(s)
    return samples

if __name__ == '__main__':
    main()

