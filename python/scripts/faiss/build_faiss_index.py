import os
from os.path import basename
import sys
import argparse
import numpy as np
import faiss

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--emb', type=str)
    parser.add_argument('--idx_dir', type=str)
    parser.add_argument('--db_samples', type=str)
    #parser.add_argument('--emb_ext', type=str, default='.emb')
    return parser.parse_args()

def main():
    args = parse_args()
    emb = args.emb
    idx_dir = args.idx_dir
    db_samples = args.db_samples
    #emb_ext = args.emb_ext

    # read database samples
    db_samples_list = read_database_samples(db_samples)

    print('Building index for: {}'.format(emb))
    # read data from file
    gt_embedding_file = emb
    embeddings_numpy = read_embeddings(gt_embedding_file, db_samples_list)

    # build faiss index for l2 distance and write to file
    l2_index = build_l2_index(embeddings_numpy)
    #base_name = '_'.join(gt_embedding.split('_')[0:2]).replace('.txt', '')
    base_name = '.'.join(basename(emb).split('.')[0:2])
    l2_index_file = idx_dir + base_name + '.idx'
    print(base_name, l2_index_file)
    faiss.write_index(l2_index, l2_index_file)

def read_database_samples(db_samples):
    # return list of database samples
    db_samples_list = []
    for line in open(db_samples, 'r'):
        db_samples_list.append(line.strip())
    return db_samples_list

def read_embeddings(gt_embedding_file, db_samples_list):
    # return array of embeddings
    embeddings = []

    # open embedding file
    for line in open(gt_embedding_file, 'r'):
        # split line
        line = line.split(' ')
        # get segment name
        sample_hap = line[0].split()[0]
        sample_ID = sample_hap.split('_')[0]
        # check if segment is in database
        if sample_ID not in db_samples_list:
            continue
        # get segment embedding
        segment_embedding = np.array([float(i) for i in line[1:]])
        # append to embeddings
        embeddings.append((sample_hap, segment_embedding))

    # convert data to numpy array and reshape
    data_array = np.array([i[1] for i in embeddings])
    data_array = data_array.reshape(data_array.shape[0], data_array.shape[1])
    return data_array

def build_l2_index(data_array):
    # build faiss index for l2 distance
    index = faiss.IndexFlatL2(data_array.shape[1])
    index.add(data_array)
    return index

if __name__ == '__main__':
    main()
