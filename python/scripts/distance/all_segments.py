import argparse
import os

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--dist_dir')
    parser.add_argument('--dist_ext')
    parser.add_argument('--num_seg')
    parser.add_argument('--hap')
    return parser.parse_args()

def main():
    args = get_args()

    samples_dict = dict()

    for f in os.listdir(args.dist_dir):
        if f.endswith(args.dist_ext):
            #print(f)
            file_path = args.dist_dir + f
            seg = int(f.split('.')[2])
            #print('adding segment ', seg)
            read_segment(samples_dict, file_path, int(args.num_seg), seg, args.hap)

    header = 'sampleID ' + ' '.join([str(i) for i in range(int(args.num_seg))])
    print(header)
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


if __name__ == '__main__':
    main()
