import argparse
import matplotlib.pyplot as plt
import sys, os

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--vcf')
    parser.add_argument('--map')
    parser.add_argument('--out_png')
    parser.add_argument('--num_vars')

    return parser.parse_args()

def main():
    args = get_args()
    bp = get_bps(args.vcf, int(args.num_vars))
    plot_hist(bp, args.out_png)

def get_bps(vcf_file, num_variants):
    variants = []
    nv = 0
    f = open(vcf_file)
    while(nv < num_variants):
        line = f.readline()
        if '#' in line:
            continue
        pos = int(line.strip().split()[1])
        variants.append(pos)
        nv += 1
    print('read ', nv, ' variants.')
    f.close()
    return variants

def plot_hist(bps, outpng):
    png_name = outpng + 'variants.png'
    plt.figure(figsize=(28, 15))
    ax1 = plt.subplot(111)
    ax1.set_title('distribution of variants', fontsize=30)
    ax1.hist(bps, bins=100, color='olivedrab')
    ax1.set_xlabel('basepair position', fontsize=20)
    ax1.spines.right.set_visible(False)
    ax1.spines.top.set_visible(False)

    plt.savefig(png_name)

if __name__ == '__main__':
    main()
