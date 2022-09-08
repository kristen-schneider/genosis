import read_vcf
import argparse

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--vcf')
    parser.add_argument('--map')
    parser.add_argument('--out_dir')

    return parser.parse_args()

def main():
    args = get_args()
    read_vcf.read_vcf_block(args.vcf, 61879, 1093148)
    

if __name__ == '__main__':
    main()
