import os
import pysvs
import argparse
import numpy as np
from os.path import basename

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--emb_file', type=str)
    parser.add_argument('--idx_dir', type=str)
    parser.add_argument('--db_samples', type=str)
    return parser.parse_args()

def main():
    args = parse_args()
    emb_file = args.emb_file
    idx_dir = args.idx_dir
    db_samples = args.db_samples

    # simple file name
    root_name = basename(emb_file)
    # read database sample IDs
    db_samples_list = read_database_samples(db_samples)
    
    # svs setup
    parameters = pysvs.VamanaBuildParameters(
        graph_max_degree = 64,
        window_size = 128,
        num_threads = 4,
    )

    # tuple (id, numpy array for embeddings)
    numpy_embeddings = read_embeddings(emb_file, db_samples_list)
    # create svs index
    index = pysvs.Vamana.build(
            parameters, 
            numpy_embeddings, 
            pysvs.DistanceType.L2,
    )
    # save svs index
    index.save(
            os.path.join(idx_dir, root_name.replace('emb', '') + 'config'),
            os.path.join(idx_dir, root_name.replace('emb', '') + 'graph'),
            os.path.join(idx_dir, root_name.replace('emb', '') + 'data'),
            )
    
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
        line = line.strip().split(' ')
        # get segment name
        sample_hap = line[0].strip()
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
