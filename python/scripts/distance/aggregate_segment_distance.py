import argparse
from collections import defaultdict
import os

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--out_dir')
    parser.add_argument('--ext')
    parser.add_argument('--num_seg')
    return parser.parse_args()

def main():
    args = get_args()
    queries_dict = defaultdict(dict)

    for f in os.listdir(args.out_dir):
        if f.endswith(args.ext):
            file_path = args.out_dir + f
            seg = int(f.split('.')[3])
            read_segment(queries_dict, file_path, int(args.num_seg), seg)
    
    for q in queries_dict:
        print('Query: ', q)
        header = 'sampleID ' + ' '.join([str(i) for i in range(int(args.num_seg))])
        print(header)
        for sample in queries_dict[q]:
            line = sample + ' ' + ' '.join([str(dist) for dist in queries_dict[q][sample]])
            print(line)



def read_segment(queries_dict, seg_dist_file, num_seg, seg):
    f = open(seg_dist_file, 'r')
    for line in f:
        if 'Query' in line:
            query = line.strip().split(':')[1].strip()
        else:
            L = line.strip().split()
            try:
                sample = L[0]
                dist = float(L[1])
            except IndexError:
                continue
            try:
                queries_dict[query][sample][seg] = dist
            except KeyError:
                try:
                    empty_list = [-1] * int(num_seg)
                    empty_list[seg] = dist
                    queries_dict[query][sample] = empty_list
                except KeyError:
                    empty_list = [-1] * int(num_seg)
                    empty_list[seg] = dist
                    queries_dict[query] = {sample: empty_list}
                #queries_dict[query] = dist_dict[sample]
            #except KeyError:
            #    dist_dict = dict()
            #    empty_list = [-1] * int(num_seg)
            #    dist_dict[sample] = empty_list
            #    dist_dict[sample][seg] = dist
            #    queries_dict[query][sample][seg] = dist
    #print(queries_dict)
    f.close()

if __name__ == '__main__':
    main()
