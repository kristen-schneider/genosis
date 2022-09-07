import argparse

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--vcf')
    parser.add_argument('--start')
    parser.add_argument('--end')
    parser.add_argument('--out_vcf')
    return parser.parse_args()

def main():
    args = get_args()
    write_vcf_segment(args.vcf, int(args.start), int(args.end), args.out_vcf)
    #for ibd_match in segment_boundaries:
    #    [ibd_start, ibd_end] = segment_boundaries[ibd_match]
    #    out_vcf_name = args.out_dir + 'ibd.' + str(ibd_match) + '.vcf'
    #    write_vcf_segment(args.vcf, ibd_start, ibd_end, out_vcf_name)

def write_vcf_segment(vcf_file, start, end, out_vcf):
    print('writing to ', out_vcf)
    o = open(out_vcf, 'w')
    f = open(vcf_file, 'r')
    for line in f:
        if '#' in line:
            o.write(line)
            continue
        L = line.strip().split()
        pos = int(L[1])
        if pos >= start and pos <= end:
            o.write(line)
    f.close()
    o.close()
    print('done writint to ', out_vcf)

if __name__ == '__main__':
    main()
