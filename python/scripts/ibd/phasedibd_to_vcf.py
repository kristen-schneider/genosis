import argparse

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--vcf')
    parser.add_argument('--ilash')
    parser.add_argument('--out_dir')
    return parser.parse_args()

def main():
    args = get_args()
    [segment_boundaries, segment_pairs, segment_length] = read_ilash_IBD(args.ilash)
    for ibd_match in segment_boundaries:
        [ibd_start, ibd_end] = segment_boundaries[ibd_match]
        out_vcf_name = args.out_dir + 'ibd.' + str(ibd_match) + '.vcf'
        write_vcf_segment(args.vcf, ibd_start, ibd_end, out_vcf_name)

    x = ''

def write_vcf_segment(vcf_file, start, end, out_vcf):
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



def read_phasedibd_IBD(ibd_file):
    segment_bp_boundaries = dict()
    segment_pairs = dict()
    segment_cm_boundaries = dict()

    ibd_i = 0

    f = open(ibd_file, 'r')
    for line in f:
        L = line.strip().split()
        sample1 = L[5]
        sample2 = L[6]
        startBP = int(L[1])
        endBP = int(L[2])
        startCM = float(L[3])
        endCM = float(L[4])

        segment_bp_boundaries[ibd_i] = [startBP, endBP]
        segment_pairs[ibd_i] = [sample1, sample2]
        segment_cm_boundaries[ibd_i] = [startCM, endCM]
        ibd_i += 1

    f.close()
    return [segment_boundaries, segment_pairs, segment_length]

if __name__ == '__main__':
    main()
