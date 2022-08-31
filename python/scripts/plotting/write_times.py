import argparse
import re
import sys, os

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--data_dir')
    parser.add_argument('--plink_ext')
    parser.add_argument('--faiss_enc_ext')
    parser.add_argument('--faiss_emb_ext')
    parser.add_argument('--out_dir')
    parser.add_argument('--num_seg')
    return parser.parse_args()

def main():
    args = get_args()

    [plink_seg_times, faiss_enc_seg_times, faiss_emb_seg_times] = read_file_times(args)
    write_file_times(plink_seg_times, args.out_dir + 'plink.seg.times', args.num_seg)
    write_file_times(faiss_enc_seg_times, args.out_dir + 'faiss_enc.seg.times', args.num_seg)
    write_file_times(faiss_emb_seg_times, args.out_dir + 'faiss_emb.seg.times', args.num_seg)


def write_file_times(times_dict, out_file, num_seg): 
    print(out_file)
    o = open(out_file, 'w')
    for s in range(int(num_seg)):
        line = str(s) + '\t' + str(times_dict[str(s)]) + '\n'
        o.write(line)
    o.close



def read_file_times(args):
    plink_seg_times = dict()
    faiss_enc_seg_times = dict()
    faiss_emb_seg_times = dict()


    # get faiss times for each segment
    for segment_file in os.listdir(args.data_dir):
        f = os.path.join(args.data_dir, segment_file)
        try: seg_i = segment_file.split('.')[2]
        except IndexError: continue


        # seg_i's plink file
        if (seg_i in f and args.plink_ext in f):
            s_time = get_plink_time(f)
            plink_seg_times[seg_i] = s_time

        # seg_i's faiss_enc file
        elif (seg_i in f and args.faiss_enc_ext in f):
            s_time = get_faiss_time(f)
            faiss_enc_seg_times[seg_i] = s_time

        # seg_i's faiss_emb file
        elif (seg_i in f and args.faiss_emb_ext in f):
            s_time = get_faiss_time(f)
            faiss_emb_seg_times[seg_i] = s_time

    return [plink_seg_times, faiss_enc_seg_times, faiss_emb_seg_times]


def get_faiss_time(faiss_file):
    # microseconds
    seg_time = -1
    f = open(faiss_file, 'r')
    for line in f:
        if re.search(r'^TIME', line):
            time_line = line.split(':')
            seg_time = int(time_line[2])
    f.close()
    return (seg_time / 1e6)

def get_plink_time(plink_file):
    # seconds
    seg_time = -1
    f = open(plink_file, 'r')
    for line in f:
        # Start time: Tue Aug 30 18:33:35 2022
        if re.search(r'^Start time:', line):
            start_minute = int(line.split()[5].split(':')[-2])
            start_second = int(line.split()[5].split(':')[-1])
        elif re.search(r'^End time:', line):
            end_minute = int(line.split()[5].split(':')[-2])
            end_second = int(line.split()[5].split(':')[-1])
    if start_minute == end_minute:
        seg_time = end_second - start_second
    elif start_minute < end_minute:
        seg_time = 60-start_second+end_second
    f.close()
    return seg_time

if __name__ == '__main__':
    main()

