import argparse
import os
import sys

data_dir = sys.argv[1]
out_dir = sys.argv[2]

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--svs_sample_results_dir')
    parser.add_argument('--sample')
    return parser.parse_args()

def main():
    args = get_args()
    svs_dir = args.svs_sample_results_dir
    all_chrm = args.sample

    sorted_dict = read_sample(svs_dir + sample)        
    topK_dict = topK(sorted_dict)
    # write to file
    with open(out_dir + sample, 'w') as f:
        for m in topK_dict:
            f.write(m + '\t' + str(topK_dict[m]) + '\n')

def read_sample(sample_file):
    hap = None
    sample_dict = {}
    with open(sample_file) as f:
        for line in f:
            if hap == None: hap = line
            elif ':' not in line: hap = line
            else:
                L = line.strip().split(':')
                match_ID = L[0].split('_')[0].strip()
                score = float(L[1].strip())
                try:
                    sample_dict[match_ID] += score
                except KeyError:
                    sample_dict[match_ID] = score       
    # sort by value
    sorted_dict = dict(sorted(sample_dict.items(), key=lambda item: item[1], reverse=True))
    return sorted_dict

def topK(sorted_dict, k=20):
    # only keep first K from sorted_dict
    topK_dict = {}
    for i, m in enumerate(sorted_dict):
        if i < k: topK_dict[m] = sorted_dict[m]
        else: break
    return topK_dict

if __name__ == '__main__':
    main()
