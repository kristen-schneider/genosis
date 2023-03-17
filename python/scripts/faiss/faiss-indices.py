import sys
import time
import numpy as np
import faiss

def main():
    segments_dir=sys.argv[1]
    encoding_file=sys.argv[2]
    base_name=sys.argv[3]
    k=int(sys.argv[4])

    # read data from file
    # print('Reading data from file: {}'.format(encoding_file))
    data = read_data(encoding_file)
    # convert data to numpy array
    data_array = np.array([i[1] for i in data])

    # L2 Index
    l2_index_file = segments_dir + base_name + '.index.l2'
    l2_out_file = segments_dir + 'results.l2'
    # build l2 index
    print('Building l2 index')
    # time index build
    start_build_time = time.time()
    l2_index = build_l2_index(data_array)
    end_build_time = time.time()
    # write index to file
    faiss.write_index(l2_index, l2_index_file)
    # search index
    print('Searching l2 index')
    # time search
    start_search_time = time.time()
    # read index from file
    l2_read_index = faiss.read_index(l2_index_file)
    D, I = l2_read_index.search(data_array, k)
    end_search_time = time.time()
    # print, and append results and timing to file
    with open(l2_out_file, 'a') as f:
        print('Segment: {}'.format(base_name), file=f)
        print('Build time: {}'.format(end_build_time - start_build_time), file=f)
        print('Search time: {}'.format(end_search_time - start_search_time), file=f)
        for i in range(len(data)):
            print('Query: {}'.format(data[i][0]), file=f)
            print('Results:', file=f)
            for j in range(k):
                print('  {}: {}'.format(data[I[i][j]][0], D[i][j]), file=f)
            print(file=f)
    # delete index
    l2_index.reset()
    l2_read_index.reset()



    # HNSW Index
    hnsw_index_file = segments_dir + base_name + '.index.hnsw'
    hnsw_out_file = segments_dir + 'results.hnsw'
    # build hnsw index
    print('Building hnsw index')
    # time index build
    start_build_time = time.time()
    hnsw_index = build_hnsw_index(data_array)
    end_build_time = time.time()
    # write index to file
    faiss.write_index(hnsw_index, hnsw_index_file)
    # search index
    print('Searching hnsw index')
    # time search
    start_search_time = time.time()
    # read index from file
    hnsw_read_index = faiss.read_index(hnsw_index_file)
    D, I = hnsw_read_index.search(data_array, k)
    end_search_time = time.time()
    # print results and timing to file
    with open(hnsw_out_file, 'a') as f:
        print('Segment: {}'.format(base_name), file=f)
        print('Build time: {}'.format(end_build_time - start_build_time), file=f)
        print('Search time: {}'.format(end_search_time - start_search_time), file=f)
        for i in range(len(data)):
            print('Query: {}'.format(data[i][0]), file=f)
            print('Results:', file=f)
            for j in range(k):
                print('  {}: {}'.format(data[I[i][j]][0], D[i][j]), file=f)
            print(file=f)
    # delete index
    hnsw_index.reset()
    hnsw_read_index.reset()


    ## IVFPQR Index
    #ivfpqr_index_file = segments_dir + base_name + '.index.ivfpqr'
    #ivfpqr_out_file = segments_dir + base_name + '.results.ivfpqr'
    ## build ivfpqr index
    #print('Building ivfpqr index')
    ## time index build
    #start_build_time = time.time()
    #ivfpqr_index = build_ivfpqr_index(data_array)
    #end_build_time = time.time()
    ## write index to file
    #faiss.write_index(ivfpqr_index, ivfpqr_index_file)
    ## search index
    #print('Searching ivfpqr index')
    ## time search
    #start_search_time = time.time()
    ## read index from file
    #ivfpqr_read_index = faiss.read_index(ivfpqr_index_file)
    #D, I = ivfpqr_read_index.search(data_array, k)
    #end_search_time = time.time()
    ## print results
    #with open(ivfpqr_out_file, 'a') as f:
    #    print('Segment: {}'.format(base_name), file=f)
    #    print('Build time: {}'.format(end_build_time - start_build_time), file=f)
    #    print('Search time: {}'.format(end_search_time - start_search_time), file=f)
    #    for i in range(len(data)):
    #        print('Query: {}'.format(data[i][0]), file=f)
    #        print('Results:', file=f)
    #        for j in range(k):
    #            print('  {}: {}'.format(data[I[i][j]][0], D[i][j]), file=f)
    #        print(file=f)
    ## delete index
    #ivfpqr_index.reset()
    #ivfpqr_read_index.reset()


def build_l2_index(data_array):
    d = data_array.shape[1]
    index = faiss.IndexFlatL2(d)
    index.add(data_array)
    return index

def build_hnsw_index(data_array):
    d = data_array.shape[1]
    index = faiss.IndexHNSWFlat(d, 32)
    index.hnsw.efConstruction = 200
    index.hnsw.efSearch = 32
    index.add(data_array)
    return index

def build_ivfpqr_index(data_array):
    d = 128
    quantizer = faiss.IndexFlatL2(d)
    n_list = 25 # number of clusters
    q = 8 # number of bits per sub-vector
    m = 8 # number of sub-vectors

    index = faiss.IndexIVFPQ(quantizer, d, n_list, m, q)
    index.train(data_array)
    index.add(data_array)
    return index


def read_data(data_file):
    data = []
    with open(data_file, 'r') as f:
        for line in f:
            # split line into id and vector
            id, vector = line.split(' ', 1)
            # convert vector to numpy array
            vector = np.fromstring(vector, sep=' ')
            # append id and vector to data
            data.append((id, vector))

    return data

if __name__ == '__main__':
    main()
