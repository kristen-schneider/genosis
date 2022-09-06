import argparse
import matplotlib.pyplot as plt
import sys
import os
import numpy as np

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--data_dir')
    parser.add_argument('--out_png')
    parser.add_argument('--num_snps')
    parser.add_argument('--chrm')

    return parser.parse_args()

def main():
    args = get_args()

    bp_distance_file = args.data_dir + str(args.num_snps) + '.basepair_distance.txt'
    bp_start_end = read_data_file(bp_distance_file)
    # cm_est_distance_file = args.data_dir + str(args.num_snps) + '.centimorgan_est_distance.txt'
    # cm_est_count = read_data_file(cm_est_distance_file)
    cm_map_distance_file = args.data_dir + str(args.num_snps) + '.centimorgan_map_distance.txt'
    cm_map_count = read_data_file(cm_map_distance_file)

    # plot_positions(bp_start_end, cm_est_count, cm_map_count, args.out_png, args.num_snps)
    plot_snp_vs_cm(bp_start_end, cm_map_count, args.out_png, args.num_snps, chrom)

def read_data_file(data_file):
    data_dict = dict()
    f = open(data_file, 'r')
    for line in f:
        L = line.strip().split(',')
        seg = int(L[0])
        dist = float(L[1])
        data_dict[seg] = dist
    f.close()
    return data_dict

def plot_snp_vs_cm(seg_bp, seg_cm_map, out_png, num_snps, chrom):
    png_name = out_png + num_snps + '.variants.png'
    plt.figure(figsize=(28, 15))
    ax1 = plt.subplot(111)
    title = 'Chr ' \
            + str(chrom) \
            + 'Lengths of Segments\n(seg size = ' + num_snps + ')'
    ax1.set_title(title, fontsize=30)
    # base pairs
    x_bp_data = []
    for seg in seg_bp:
        x_bp_data.append(seg_bp[seg])
    # cm pairs map
    y_cm_map_data = []
    for seg in seg_cm_map:
        y_cm_map_data.append(seg_cm_map[seg])

    # linear regression
    coef = np.polyfit(x_bp_data, y_cm_map_data, 1)
    poly1d_fn = np.poly1d(coef)

    bp = ax1.plot(x_bp_data, y_cm_map_data, 'yo',
                     x_bp_data, poly1d_fn(x_bp_data), '--k')
    # bp = ax1.scatter(x_bp_data, y_cm_map_data, color='olivedrab')
    ax1.spines.top.set_visible(False)
    ax1.spines.right.set_visible(False)
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
