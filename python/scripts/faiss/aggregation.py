import argparse
from collections import defaultdict
import os

import read_faiss
import utils

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--out_dir')
    parser.add_argument('--ext')
    parser.add_argument('--num_segments')
    parser.add_argument('--query')
    return parser.parse_args()

def main():
    args = get_args()
    
    i = 0
    # get distances data structure
    queries_dict = defaultdict(dict)
    for f in os.listdir(args.out_dir):
        if f.endswith(args.ext):
            i+=1
            file_path = args.out_dir + f
            seg = int(f.split('.')[3])
            print(f)
            read_faiss.make_query_dict(queries_dict, file_path, int(args.num_segments), seg)
            
            #if i == 4: break

    # write distances
    for q in queries_dict:
        o = open(q+'.txt', 'w')
        o.write('sampleID, dist (by segment)\n')
        for sample_id in queries_dict[q]:
            #sum_distances = utils.sum_distances(queries_dict[q][sample_id])
            line = sample_id + ' ' + ' '.join([d for d in queries_dict[q][sample_id]]) + '\n'
            o.write(line)
        o.close()
    # sum distances
    #for q in queries_dict:
    #    for sample_id in queries_dict[q]:
    #        sum_distances = utils.sum_distances(queries_dict[q][sample_id])
    #        print(sample_id, sum_distances)
            

if __name__ == '__main__':
    main()
