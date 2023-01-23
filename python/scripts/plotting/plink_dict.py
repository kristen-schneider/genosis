import argparse
from collections import defaultdict

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--plink_file')
    parser.add_argument('--out_file')
    return parser.parse_args()

def main():
    args = get_args()
    plink_file = args.plink_file
    out_file = args.out_file

    plink_dict = get_plink_dict(plink_file)
    write_plink_dict(plink_dict, out_file)

def write_plink_dict(plink_dict, out_file):
    o = open(out_file, 'w')
    o.write('sample1_ID sample2_ID dist\n')
    for sample1 in plink_dict:
        for sample2 in plink_dict[sample1]:
            line = sample1 + ' ' + sample2 + ' ' + str(plink_dict[sample1][sample2]) + '\n'
            o.write(line)
    o.close()

def get_plink_dict(plink_file):
    plink_dict = defaultdict(dict)

    f = open(plink_file, 'r')
    header = f.readline()
    for line in f:
        l = line.strip().split()
        ID1 = l[1]
        ID2 = l[3]
        try:
            plink_distance = float(l[11])
        except ValueError:
            plink_distance = -1
        plink_dict[ID1][ID2] = plink_distance
        plink_dict[ID2][ID1] = plink_distance
    f.close()
    return plink_dict

if __name__ == '__main__':
    main()
