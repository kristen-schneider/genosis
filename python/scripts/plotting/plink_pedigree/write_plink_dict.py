import argparse
import sys

sys.path.insert(0, '../..')
from ibd import read_plink

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--plink')
    parser.add_argument('--out')
    return parser.parse_args()

def main():
    args = get_args()
    print('Reading Plink dict...')
    plink_dict = read_plink.make_plink_pairs_dict(args.plink)
    print('Writing Plink dict...')
    write_plink_dict(plink_dict, args.out)

def write_plink_dict(plink_dict, out_file):
    o = open(out_file, 'w')
    o.write('sample1_ID sample2_ID dist\n')
    for sample1 in plink_dict:
        for sample2 in plink_dict[sample1]:
            line = sample1 + ' ' + sample2 + ' ' + str(plink_dict[sample1][sample2]) + '\n'
            o.write(line)
    o.close()

if __name__ == '__main__':
    main()
