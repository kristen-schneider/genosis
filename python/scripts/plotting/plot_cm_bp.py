import argparse
import matplotlib.pyplot as plt
import sys
import os
import numpy as np

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--bp_endpoints')
    parser.add_argument('--cm_distances')
    parser.add_argument('--out_png')
    parser.add_argument('--chrm')

    return parser.parse_args()

def main():
    args = get_args()

    snp_sizes = [1000, 5000, 10000, 50000, 100000, 1000000]
    bp_files = [args.bp_endpoints + str(s) for s in snp_sizes]
    cm_files = [args.cm_distances + str(s) for s in snp_sizes]
    for size in range(len(snp_sizes)):
        basepair_dict = read_bp_file(bp_files[size])
        cm_dict = read_cm_file(cm_files[size])
        plot_bp_vs_cm(basepair_dict, cm_dict, args.out_png, snp_sizes[size], args.chrm)

def read_cm_file(data_file, delim=' '):
    data_dict = dict()
    f = open(data_file, 'r')
    header = f.readline()
    for line in f:
        L = line.strip().split(delim)
        seg = int(L[0])
        dist = float(L[1])
        data_dict[seg] = dist
    f.close()
    return data_dict

def read_bp_file(data_file, delim=' '):
    data_dict = dict()
    f = open(data_file, 'r')
    header = f.readline()
    for line in f:
        L = line.strip().split(delim)
        seg = int(L[0])
        bp_start = int(L[1])
        bp_end = int(L[2])
        data_dict[seg] = bp_end - bp_start
    f.close()
    return data_dict

def plot_bp_vs_cm(basepair_dict, cm_dict, out_png, num_snps, chrom):
    png_name = out_png + 'chr' + str(chrom) + '.' \
               + str(num_snps) + '.cM-bps.png'
    title = 'Relationship between base pairs and centimorgans' \
            '\n(Chromosome ' \
            + str(chrom) \
            + ', segment size = ' \
            + str(num_snps) + ' SNPs)'

    # base pairs
    x_bp_data = []
    for seg in basepair_dict:
        x_bp_data.append(basepair_dict[seg])
    # centimorgans
    y_cm_map_data = []
    for seg in cm_dict:
        y_cm_map_data.append(cm_dict[seg])


    # linear regression
    b, a = np.polyfit(x_bp_data, y_cm_map_data, deg=1)
    xseq = np.linspace(0, len(x_bp_data), num=len(x_bp_data))


    plt.figure(figsize=(15, 12))
    ax1 = plt.subplot(111)
    ax1.set_title(title, fontsize=30)
    ax1.spines.top.set_visible(False)
    ax1.spines.right.set_visible(False)
    bp = ax1.plot(x_bp_data, y_cm_map_data, 'yo',
                  x_bp_data, a + b * xseq, ':k')
    ax1.set_xscale('log')
    ax1.set_xlabel('basepairs', fontsize=20)
    ax1.set_ylabel('centimorgans', fontsize=20)

    plt.savefig(png_name)

def plot_positions(seg_bp, seg_cm_estimate, seg_cm_map, out_png, num_snps):
    png_name = out_png + num_snps + '.variants.png'
    plt.figure(figsize=(28, 15))
    ax1 = plt.subplot(111)
    ax1.set_title('Lengths of Segments', fontsize=30)

    x = range(len(seg_bp))
    # base pairs
    y_bp_data = []
    for seg in seg_bp:
        y_bp_data.append(seg_bp[seg])
    # cm pairs estimate
    y_cm_estimate_data = []
    for seg in seg_cm_estimate:
        y_cm_estimate_data.append(seg_cm_estimate[seg])
    # cm pairs map
    y_cm_map_data = []
    for seg in seg_cm_map:
        y_cm_map_data.append(seg_cm_map[seg])

    ax1.set_xticks(x)
    bp = ax1.plot(x, y_bp_data, color='olivedrab', label='Base Pairs')
    ax1.spines.top.set_visible(False)
    ax1.spines.right.set_visible(False)
    ax1.set_xlabel('Segment', fontsize=20)
    ax1.set_ylabel('Length (bp)', fontsize=20)

    ax2 = ax1.twinx()
    cm_e = ax2.plot(x, y_cm_estimate_data, color='firebrick', label='Centimorgan Estimated')
    ax2.set_ylabel('Length (cM)', fontsize=20)
    ax2.spines.top.set_visible(False)
    ax2.spines.left.set_visible(False)
    cm_m = ax2.plot(x, y_cm_map_data, color='steelblue', label='Centimorgan Mapped')
    ax2.spines.left.set_visible(False)
    ax2.spines.top.set_visible(False)

    lns = bp + cm_e + cm_m
    labs = [l.get_label() for l in lns]
    ax2.legend(lns, labs, fontsize=40)

    plt.savefig(png_name)

if __name__ == '__main__':
    main()
