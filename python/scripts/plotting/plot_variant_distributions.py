import argparse
import matplotlib.pyplot as plt
import sys, os

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--map')
    parser.add_argument('--out_dir')
    parser.add_argument('--cm')
    parser.add_argument('--chr')
    return parser.parse_args()

def main():
    args = get_args()
    cm_boundaries = get_cm_boundaries(args.map, int(args.cm))
    write_cm_boundaries(cm_boundaries, args.chr, args.cm, args.out_dir)



    # [variant_distribution, cm_length] = get_variant_counts(args.map, int(args.cm))
    # plot_hist(variant_distribution, args.cm, args.chr, cm_length, args.out_png)

def write_cm_boundaries(cm_boundaries, chrm, cm, out_dir):
    out_file = out_dir + 'chr' + str(chrm) + '.' + str(cm) + 'cm-segment-boundaries'
    f = open(out_file, 'w')
    f.write('segment, start_bp, end_bp\n')
    for seg in cm_boundaries:
        start = cm_boundaries[seg][0]
        end = cm_boundaries[seg][1]
        line = str(seg) + ',' + str(start) + ',' + str(end) + '\n'
        f.write(line)
    f.close()

def get_cm_boundaries(map_file, cm_length):
    cm_boundaries = dict()
    v_cm = 0
    cm_i = 0
    f = open(map_file, 'r')
    seg_bps = []
    cm_i_length = (cm_i + cm_length)

    while (v_cm < cm_i_length):
        L = f.readline().strip().split()
        try:
            chr = int(L[0])
            v_cm = float(L[2])
            v_pos = int(L[3])
            seg_bps.append(v_pos)
        except IndexError:
            if len(L) == 0:
                break
        if v_cm >= cm_i_length:
            cm_boundaries[cm_i] = [seg_bps[0], seg_bps[-1]]
            seg_bps = []
            cm_i += 1
            cm_i_length = cm_i + cm_length

    # last cm segment
    cm_boundaries[cm_i] = [seg_bps[0], seg_bps[-1]]
    f.close()
    return cm_boundaries

def get_variant_counts(map_file, cm_length):
    variants = dict()
    variant_count = 0
    v_cm = 0
    cm_i = 0
    f = open(map_file, 'r')
    cm_i_length = (cm_i + cm_length)

    while(v_cm < cm_i_length):
        L = f.readline().strip().split()
        try:
            chr = int(L[0])
            v_cm = float(L[2])
            v_pos = int(L[3])
            variant_count += 1
        except IndexError:
            if len(L) == 0:
                break
        if v_cm >= cm_i_length:
            variants[cm_i] = variant_count
            cm_i += 1
            cm_i_length = cm_i + cm_length

    # last cm segment
    variants[cm_i] = variant_count
    cm_i += 1

    cm_total = v_cm
    f.close()
    return [variants, cm_total]

def plot_hist(variant_distribution, cm, chr, cm_length, out_png):
    png_name = out_png + 'chr' + chr + '.variants_distribution.png'
    png_title = 'Distribution of variants for ' \
                + str(cm) \
                + ' cM segments' \
                + '\n(Chromosome ' + str(chr) + ')' \
                + '\n Full Distance = ' + str(cm_length) + 'cM'

    plot_data = []
    for cm in variant_distribution:
        plot_data.append(variant_distribution[cm])

    plt.figure(figsize=(28, 15))
    ax1 = plt.subplot(111)
    ax1.set_title(png_title, fontsize=30)
    ax1.hist(plot_data, bins=30, color='olivedrab')
    ax1.set_xlabel('number of SNPs', fontsize=20)
    ax1.spines.right.set_visible(False)
    ax1.spines.top.set_visible(False)

    plt.savefig(png_name)

if __name__ == '__main__':
    main()
