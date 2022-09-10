import argparse
import matplotlib.pyplot as plt
import numpy as np

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--cm_snps')
    parser.add_argument('--out_png')
    parser.add_argument('--chrm')

    return parser.parse_args()

def main():
    args = get_args()

    cm_sizes = [1, 2, 3]
    all_dict_data = dict()

    xrange = [0, 300000]
    cm_files = [args.cm_snps + str(s) for s in cm_sizes]
    for size in range(len(cm_sizes)):
        seg_dict = read_cm_file(cm_files[size])
        plot_single_snp_vs_cm(seg_dict, args.out_png, cm_sizes[size], args.chrm, xrange)
        all_dict_data[size] = seg_dict

    plot_all_snp_vs_cm(all_dict_data, args.out_png, args.chrm, xrange)

def read_cm_file(data_file, delim=' '):
    data_dict = dict()
    f = open(data_file, 'r')
    header = f.readline()
    for line in f:
        L = line.strip().split(delim)
        seg = int(L[0])
        SNPs = int(L[1])
        data_dict[seg] = SNPs
    f.close()
    return data_dict

def plot_single_snp_vs_cm(seg_dict, out_png, cm_size, chrom, xrange):
    png_name = out_png + 'chr' + str(chrom) + '.' \
               + str(cm_size) + '.cM-SNPs.png'

    title = 'SNP counts for ' \
            '\n(Chromosome ' \
            + str(chrom) \
            + ', segment size = ' \
            + str(cm_size) + ' cMs)'
    snp_data = []
    for seg in seg_dict:
        snp_data.append(seg_dict[seg])

    plt.figure(figsize=(15, 12))
    ax1 = plt.subplot(111)
    ax1.set_title(title, fontsize=30)
    ax1.hist(snp_data, color='salmon')
    ax1.spines.top.set_visible(False)
    ax1.spines.right.set_visible(False)
    ax1.set_xlabel('SNPs in segment', fontsize=20)
    ax1.set_xlim(xrange)
    plt.savefig(png_name)

def plot_all_snp_vs_cm(all_dict_data, out_png, chrom, xrange):
    png_name = out_png + 'chr' + str(chrom) \
               + '.all.cM-SNPs.png'

    all_hist_data = []
    for d in all_dict_data:
        snp_data = []
        for seg in all_dict_data[d]:
            snp_data.append(all_dict_data[d][seg])
        all_hist_data.append((snp_data))

    plt.figure(figsize=(15, 12))

    ax1 = plt.subplot(311)
    ax1.set_title('1 cenitmorgan', loc='right', fontsize=25)
    ax1.hist(all_hist_data[0], color='salmon')
    ax1.set_xlim(xrange)
    ax1.set_xticks([])
    ax1.spines.top.set_visible(False)
    ax1.spines.right.set_visible(False)
    ax1.spines.bottom.set_visible(False)

    ax2 = plt.subplot(312)
    ax2.set_title('2 cenitmorgans', loc='right', fontsize=25)
    ax2.hist(all_hist_data[1], color='salmon')
    ax2.set_xlim(xrange)
    ax2.set_xticks([])
    ax2.spines.top.set_visible(False)
    ax2.spines.right.set_visible(False)
    ax2.spines.bottom.set_visible(False)

    ax3 = plt.subplot(313)
    ax3.set_title('3 cenitmorgans', loc='right', fontsize=25)
    ax3.hist(all_hist_data[2], color='salmon')
    ax3.set_xlim(xrange)
    ax3.spines.top.set_visible(False)
    ax3.spines.right.set_visible(False)
    ax3.set_xlabel('SNPs in segment', fontsize=20)


    plt.savefig(png_name)

if __name__ == '__main__':
    main()
