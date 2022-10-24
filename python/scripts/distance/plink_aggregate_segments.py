import argparse
from collections import defaultdict
import os

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--plink')
    #parser.add_argument('--dist_dir')
    #parser.add_argument('--dist_ext')
    #parser.add_argument('--num_seg')
    parser.add_argument('--query')
    return parser.parse_args()

def main():
    args = get_args()

    samples_dict = dict()
    plink_dict = make_plink_pairs_dict(args.plink)
    query_dict = plink_dict[args.query]
    header = 'sampleID dist'#' '.join([str(i) for i in range(int(args.num_seg))])
    print(header)
    for sample in query_dict:
        line = sample + ' ' + str(query_dict[sample])#' '.join([str(dist) for dist in samples_dict[sample]])
        print(line)

def agg_plink():
    stop = 0
    for f in os.listdir(args.dist_dir):
        if f.endswith(args.dist_ext):
            file_path = args.dist_dir + f
            seg = int(f.split('.')[3])
            plink_seg_dict = make_plink_pairs_dict(args.dist_dir+f)
            plink_query = args.query.split('_')[1]
            query_plink = plink_seg_dict[plink_query]
            for sample in query_plink:
                try:
                    samples_dict[sample][seg] = query_plink[sample]
                except KeyError:
                    samples_dict[sample] = [-1] * int(args.num_seg)
                    samples_dict[sample][seg] = query_plink[sample]
            #samples_dict[seg] = plink_seg_dict[args.query]
            #read_segment(samples_dict, file_path, int(args.num_seg), seg, args.hap)
            stop += 1
        #if stop == 1:
        #    break

    header = 'sampleID ' + ' '.join([str(i) for i in range(int(args.num_seg))])
    print(header)
    for sample in samples_dict:
        line = sample + ' ' + ' '.join([str(dist) for dist in samples_dict[sample]])
        print(line)

def make_plink_pairs_dict(plink_file, delim=' '):
    '''
    reads output from plink genome into a dictionary
    whos key is a sample name and whose value is another
    dictionary. the second dictionary has key as sample
    match (pair) and has value as the plink dist

    '''
    plink_pairs_dict = defaultdict(dict)

    f = open(plink_file, 'r')
    header = f.readline()
    for line in f:
        L = line.strip().split()
        sample1 = L[1]
        sample2 = L[3]
        dist = float(L[11])
        plink_pairs_dict[sample1][sample2] = dist
        plink_pairs_dict[sample2][sample1] = dist
    f.close()
    return plink_pairs_dict

if __name__ == '__main__':
    main()
