import argparse
import os

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--query_file')
    parser.add_argument('--out_dir')
    parser.add_argument('--dist_ext')
    parser.add_argument('--num_seg')
    parser.add_argument('--hap')
    return parser.parse_args()

def main():
    args = get_args()

    samples_dict = dict()

    for f in os.listdir(args.out_dir):
        if f.endswith(args.dist_ext):
            file_path = args.out_dir + f
            seg = int(f.split('.')[3])
            #print('adding segment ', seg)
            read_segment(samples_dict, file_path, int(args.num_seg), seg, args.hap)

    header = 'sampleID ' + ' '.join([str(i) for i in range(int(args.num_seg))])
    print(header)
    #print(samples_dict)
    for sample in samples_dict:
        line = sample + ' ' + ' '.join([str(dist) for dist in samples_dict[sample]])
        print(line)

def read_segment(samples_dict, seg_dist_file, num_seg, seg, hap):
    f = open(seg_dist_file, 'r')
    query_header = f.readline()
    header = f.readline()
    for line in f:
        L = line.strip().split()
        sample = L[0]
        dist = float(L[1])
        try:
            samples_dict[sample][seg] = dist
        except KeyError:
            empty_list = [-1] * int(num_seg)
            samples_dict[sample] = empty_list
            samples_dict[sample][seg] = dist
        #if '_'+hap in sample:
        #    try:
        #        samples_dict[sample][seg] = dist
        #    except KeyError:
        #        empty_list = [-1] * int(num_seg)
        #        samples_dict[sample] = empty_list
        #        samples_dict[sample][seg] = dist
    f.close()
    #if seg == 0:
    #    print(samples_dict['HG00403_0'])

def read_queries(queries_file):
    '''
    read a file of sample IDs for queries
    and return a list
    '''
    queries = []
    f = open(queries_file, 'r')
    for line in f:
        q = line.strip()
        queries.append(q)
    f.close()
    return queries

if __name__ == '__main__':
    main()
