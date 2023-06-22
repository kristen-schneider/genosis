# [imports]
import os
import pysvs
import argparse
import numpy as np
# [imports]

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--emb_dir', type=str)
    parser.add_argument('--idx_dir', type=str)
    parser.add_argument('--db_samples', type=str)
    parser.add_argument('--emb_ext', type=str, default='.emb')
    return parser.parse_args()

def main():
    args = parse_args()
    emb_dir = args.emb_dir
    idx_dir = args.idx_dir
    db_samples = args.db_samples
    emb_ext = args.emb_ext

    # read database sample IDs
    db_samples_list = read_database_samples(db_samples)

    parameters = pysvs.VamanaBuildParameters(
        graph_max_degree = 64,
        window_size = 128,
        num_threads = 4,
    )

    # for a directory of genotype embeddings, build index for l2 distance
    # only build index for database samples
    for gt_embedding in os.listdir(emb_dir):
        # check if file is a gt embedding
        if gt_embedding.endswith(emb_ext):
            print('Building index for: {}'.format(gt_embedding))
            

            # read data from file
            gt_embedding_file = emb_dir + gt_embedding
            # tuple (id, numpy array for embeddings)
            embeddings_numpy = read_embeddings(gt_embedding_file, db_samples_list)

            index = pysvs.Vamana.build(
                    parameters, 
                    embeddings_numpy, 
                    pysvs.DistanceType.L2,
            )
            print('done building index')

            index.save(
                    os.path.join(idx_dir, gt_embedding.replace('.gt', '') + '_config'),
                    os.path.join(idx_dir, gt_embedding.replace('.gt', '') + '_graph'),
                    os.path.join(idx_dir, gt_embedding.replace('.gt', '') + '_data'),
                    )
            
            ## build faiss index for l2 distance and write to file
            #l2_index = build_l2_index(embeddings_numpy)
            #base_name = '_'.join(gt_embedding.split('_')[0:2]).replace('.txt', '')
            #l2_index_file = idx_dir + base_name + '.index.l2'
            #faiss.write_index(l2_index, l2_index_file)
        # break

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
        sample_hap = line[0].split()[0]
        sample_ID = sample_hap.split('_')[0]
        # check if segment is in database
        if sample_ID not in db_samples_list:
            continue
        else:
            # get segment embedding
            segment_embedding = [float(i) for i in line[1:]]
            # append to embeddings
            embeddings.append((sample_hap, segment_embedding))
            
    # convert data to numpy array and reshape
    #np.vstack([data_array, embeddings[i][1]])
    data_array = [np.array(i[1], dtype=np.float32) for i in embeddings]
    #data_array = np.array([i[1] for i in embeddings])
    #data_array = data_array.reshape(data_array.shape[0], data_array.shape[1])
    return data_array

if __name__ == '__main__':
    main()
