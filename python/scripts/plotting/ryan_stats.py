import argparse
import re
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--plink')
    parser.add_argument('--faiss')
    parser.add_argument('--ids')
    parser.add_argument('--outfile')
    return parser.parse_args()

def main():
    args = get_args()

    idx_to_id = {}
    sample_ids = []
    with open(args.ids) as f:
        for l in f:
            idx_to_id[len(idx_to_id)] = l.rstrip()
            sample_ids.append(l.rstrip())

    query_id = None
    idxs = None
    dists = None
    Q = {}
    with open(args.faiss) as f:
        for l in f:
            if len(l) == 1:
                continue
            elif re.search(r'^TIME', l):
                continue
            elif re.search(r'^QUERY', l):
                query_id = int(l.rstrip().split()[1])
            elif query_id is not None:
                A = l.rstrip().split()
                idx = int(A[0])
                dist = float(A[1])
                if idx_to_id[query_id] not in Q:
                    Q[idx_to_id[query_id]] = []
                Q[idx_to_id[query_id]].append((idx_to_id[idx],dist))

    header_len = 1
    P = {}
    with open(args.plink) as f:
        for l in f:
            if header_len > 0:
                header_len = header_len - 1
                continue
            A = l.rstrip().split()
            a = A[1]
            b = A[3]
            dist = float(A[11])
            if a not in P: P[a] = []
            P[a].append((b,dist))

            if b not in P: P[b] = []
            P[b].append((a,dist))

    for p in P:
        P[p].sort(key=lambda tup: tup[1], reverse=True)


    i = 0
    B = []
    TPS = []
    FPS = []
    FNS = []
    for s in sample_ids:
        p_ids = [x[0] for x in P[s]]

        b = []

        for q in Q[s][:20]:
            if q[0] == s: continue
            idx = p_ids.index(q[0])
            #print(idx, q, P[s][idx])
            b.append(idx)
        i+=1
        B.append(b)
        tps = len( [x for x in b if x <= 20 ] ) 
        print(tps, b)

        TPS.append( tps )


        if i == 50: break

    fig, ax = plt.subplots()
    fig.set_size_inches(10,4)
    ax.boxplot(B)
    ax.set_ylabel('Plink rank')
    ax.set_xlabel('Top 20 FAISS hits')
    ax.spines['bottom'].set_visible(False)
    ax.set_xticklabels([])
    ax.set_xticks([])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    fig.savefig(args.outfile)

    fig, ax = plt.subplots()
    fig.set_size_inches(10,4)
    ax.hist(TPS)
    ax.set_ylabel('True positives')
    fig.savefig('tp-' + args.outfile)




if __name__ == '__main__': main()
