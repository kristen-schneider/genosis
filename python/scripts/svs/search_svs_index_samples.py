import os
import pysvs
import argparse
import numpy as np
import os.path
from collections import defaultdict
from os.path import basename

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--idx_dir', type=str, help='directory of all embedding files')
    parser.add_argument('-e', '--emb_dir', type=str, help='directory of all embedding files')
    parser.add_argument('-d', '--db_samples', type=str, help='file of all samples in database (haplotypes)')
    parser.add_argument('-q', '--query_samples', type=str, help='file of all query samples (haplotypes)')
    parser.add_argument('-k', '--knn', type=int, default=20, help='k used for knn')
    parser.add_argument('-o', '--out_dir', type=str, help='directory to write svs search results')
    return parser.parse_args()

def main():
    # get arguments from argparse
    args = parse_args()
    idx_dir = args.idx_dir
    emb_dir = args.emb_dir
    db_samples = args.db_samples
    query_samples = args.query_samples
    k = args.knn
    out_dir = args.out_dir

    # segment top hits for all queries
    # dict: key = query, value = dict of matches
    segment_hits = defaultdict(dict)

    # read database and query samples into lists
    database_samples_list = read_samples(db_samples)
    query_samples_list = read_samples(query_samples)

    # iterate through all segments in index directory (for all queries)
    print('getting results...')
    i = 0
    for seg_idx in os.listdir(idx_dir):
        if seg_idx.endswith('.config'):
            print(i+1, seg_idx)
            i += 1
            segment_hits = single_segment_knn(idx_dir+seg_idx,
                                            segment_hits,
                                            emb_dir,
                                            k,
                                            query_samples_list,
                                            database_samples_list,
                                            out_dir)

    # sort results by number of hits for each query
    for query in segment_hits:
        segment_hits[query] = {k: v for k, v in sorted(segment_hits[query].items(), key=lambda item: item[1], reverse=True)}
        

    # write out results
    print('writing output...')
    for query in segment_hits:
        print(query)
        results_file = out_dir + query + '.knn'
        rf = open(results_file, 'w')
        rf.write("QUERY: " + query + '\n')
        i = 0;
        for match in segment_hits[query]:
            if i < 20:
                #rf.write(str(database_samples_list[match]) + '\t' + str(segment_hits[query][match]) + '\n')
                rf.write(str(match) + '\t' + str(segment_hits[query][match]) + '\n')
            else: break
            i += 1
    rf.close()

def single_segment_knn(seg_idx,
                       segment_hits,
                       emb_dir,
                       k,
                       query_samples_list,
                       database_samples_list,
                       out_dir):
    # simple file name
    base = seg_idx.split('/')[-1].replace('.config', '')
    chrm, segment = base.split('.')
    root_name = base

    # get embeddings for segment
    embedding_file = emb_dir + chrm + '.' + segment + '.emb'
    # db_embeddings = read_embeddings(embedding_file, database_samples_list)
    q_embeddings = read_embeddings(embedding_file, query_samples_list)

    # get svs index from config file
    index = pysvs.Vamana(
            seg_idx,
            pysvs.GraphLoader(seg_idx.replace('.config', '.graph')),
            pysvs.VectorDataLoader(seg_idx.replace('.config', '.data'), pysvs.DataType.float32),
            pysvs.DistanceType.L2,
            num_threads = 4)

    # search index
    # I : index
    # D : distance matrix
    I, D = index.search(q_embeddings, k)

    # for all queries in our list, fill out segment_hits
    for query_idx in range(len(I)):
        #print("QUERY: ", query_samples_list[query_idx])
        q = I[query_idx]
        for match_idx in range(len(q)):
            match_ID = database_samples_list[q[match_idx]]
            try:
                #segment_hits[query_samples_list[query_idx].split('_')[0]][q[match_idx]] += 1
                segment_hits[query_samples_list[query_idx].split('_')[0]][match_ID.split('_')[0]] += 1
            except KeyError:
                try:
                    #segment_hits[query_samples_list[query_idx].split('_')[0]][q[match_idx]] = 1
                    segment_hits[query_samples_list[query_idx].split('_')[0]][match_ID.split('_')[0]] = 1
                except KeyError:
                    #segment_hits[query_samples_list[query_idx].split('_')[0]].update({q[match_idx]: 1})
                    segment_hits[query_samples_list[query_idx].split('_')[0]].update({match_ID.split('_')[0]: 1})

    return segment_hits

def read_samples(samples):
    """
    return list of samples
    """
    samples_list = []
    for line in open(samples, 'r'):
        samples_list.append(line.strip())
    return samples_list

def read_embeddings(gt_embedding_file, db_samples_list):
    """
    return array of embeddings for all dabatase samples
    """
    embeddings = []

    # open embedding file
    for line in open(gt_embedding_file, 'r'):
        # split line
        line = line.strip().split(' ')
        # get segment name
        sample_hap = line[0].split()[0]
        sample_ID = sample_hap.split('_')[0]
        # check if segment is in database
        if sample_hap not in db_samples_list:
            continue
        else:
            # get segment embedding
            segment_embedding = [float(i) for i in line[1:]]
            # append to embeddings
            embeddings.append((sample_hap, segment_embedding))
    # convert data to numpy array and reshape
    data_array = [np.array(i[1], dtype=np.float32) for i in embeddings]
    return data_array

if __name__ == '__main__':
    main()
