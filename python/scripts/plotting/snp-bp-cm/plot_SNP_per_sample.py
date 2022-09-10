import argparse
import matplotlib.pyplot as plt
import sys, os

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--base_name')
    parser.add_argument('--out_png')
    parser.add_argument('--cm')
    parser.add_argument('--chr')
    return parser.parse_args()

def main():
    args = get_args()
    all_segments = []
    segs = [12, 26, 69, 117]
    bins = range(0, 12000, 200)
    xrange = [0, 12000]
    files = [args.base_name + str(s) for s in segs]

    for f in range(len(files)):
        snp_counts = read_data_file(files[f])
        plot_single_segment(snp_counts, segs[f], args.chr, args.out_png, bins, xrange)
        all_segments.append(snp_counts)

    plot_separate_segments(all_segments, args.chr, args.cm, args.out_png, bins, xrange)

    # for s in range(0, num_segs):
    #     print('segment ', s)
    #     sample_snps_dict = count_sample_variants(cm_boundaries, s, args.vcf)
    #     segment_counts = list(sample_snps_dict.values())
    #     max_snps[s] = max(segment_counts)
    #     plot_hist(sample_snps_dict, s, args.chr, args.cm, args.out_png)
        # print('counting samples for segment ', s)
        # segment_snps = read_vcf(cm_boundaries, s, args.vcf)
        # plot_hist(segment_snps, s, args.chr, args.cm, args.out_png)

    # for seg in max_snps:
    #     print (seg, ', ', max_snps[seg])

def read_data_file(data_file, delim=' '):
    data_dict = dict()
    f = open(data_file)
    header = None
    for line in f:
        if header == None:
            header = line
            continue
        L = line.strip().split(delim)
        sample = L[0]
        SNPs = int(L[1])
        data_dict[sample] = SNPs
    f.close()
    return data_dict

def plot_single_segment(single_segment, seg, chrm, out_png, bins, xrange):
    png_name = out_png + 'chr' + chrm \
               + '.seg.' + str(seg) \
               + '.variants_distribution.png'
    png_title = 'Distribution of SNP counts per sample' \
                + '\n(Chromosome ' + str(chrm) \
                + ', segment ' + str(seg) + ')'
    single_data = []
    for sample in single_segment:
            single_data.append(single_segment[sample])

    plt.figure()#figsize=(20, 10))
    ax1 = plt.subplot(111)
    ax1.set_title(png_title)
    ax1.hist(single_data, bins=bins, color='olivedrab')
    ax1.set_xlim(xrange)
    # ax1.set_ylim(yrange)
    ax1.spines.right.set_visible(False)
    ax1.spines.top.set_visible(False)
    plt.savefig(png_name)

def plot_separate_segments(all_segments, chrm, cm, out_png, bins, xrange):
    png_name = out_png + 'chr' + chrm \
               + '.seg.all.variants_distribution.png'
    png_title = 'Distribution of SNP counts per sample' \
                + '\n(Chromosome ' + str(chrm) \
                + ', Segment Distance = ' + str(cm) + 'cM'

    all_segment_data = []
    for segment in all_segments:
        segment_hist_data = []
        for sample in segment:
            segment_hist_data.append(segment[sample])
        all_segment_data.append(segment_hist_data)

    plt.figure(figsize=(10, 15))

    ax1 = plt.subplot(411)
    ax1.set_title('Segment 12', loc='right')
    ax1.hist(all_segment_data[0], bins=bins, color='olivedrab')
    ax1.set_xlim(xrange)
    # ax1.set_ylim(yrange)
    ax1.set_xticks([])
    ax1.spines.right.set_visible(False)
    ax1.spines.top.set_visible(False)
    ax1.spines.bottom.set_visible(False)

    ax2 = plt.subplot(412)
    ax2.set_title('Segment 26', loc='right')
    ax2.hist(all_segment_data[1], bins=bins, color='olivedrab')
    ax2.set_xlim(xrange)
    # ax2.set_ylim(yrange)
    ax2.set_xticks([])
    ax2.spines.right.set_visible(False)
    ax2.spines.top.set_visible(False)
    ax2.spines.bottom.set_visible(False)

    ax3 = plt.subplot(413)
    ax3.set_title('Segment 69', loc='right')
    ax3.hist(all_segment_data[2], bins=bins, color='olivedrab')
    ax3.set_xlim(xrange)
    # ax3.set_ylim(yrange)
    ax3.set_xticks([])
    ax3.spines.right.set_visible(False)
    ax3.spines.top.set_visible(False)
    ax3.spines.bottom.set_visible(False)

    ax4 = plt.subplot(414)
    ax4.set_title('Segment 117', loc='right')
    ax4.hist(all_segment_data[3], bins=bins, color='olivedrab')
    ax4.set_xlabel('number of SNPs per sample', fontsize=15)
    ax4.set_xlim(xrange)
    # ax4.set_ylim(yrange)
    ax4.spines.right.set_visible(False)
    ax4.spines.top.set_visible(False)

    # plt.subplots_adjust(bottom=0.1,
    #                     top=0.9,
    #                     wspace=0.0,
    #                     hspace=0.3)
    plt.savefig(png_name)



#-------------------------------------------------#



def read_vcf(cm_boundaries, seg_q, vcf_file):

    query_boundaries = cm_boundaries[seg_q]
    start = query_boundaries[0]
    end = query_boundaries[1]

    segment_snps = []
    f = open(vcf_file, 'r')
    for line in f:
        # header
        if '#' in line:
            continue
        L = line.strip().split()
        pos = int(L[1])
        if (pos >= start and pos <= end):
            genotypes = L[9:]
            num_samples = count_samples(genotypes)
            segment_snps.append(num_samples)
    f.close()
    return segment_snps

def count_samples(genotypes):
    variant_samples = 0
    for gt in genotypes:
        if '1' in gt or '2' in gt:
            variant_samples += 1
    return variant_samples

def read_boundary_data(boundary_file, delim=' '):
    boundary_dict = dict()
    f = open(boundary_file)

    header = None
    for line in f:
        if header == None:
            header = line
            continue
        L = line.strip().split(delim)
        seg = int(L[0])
        start = int(L[1])
        end = int(L[2])
        boundary_dict[seg] = [start, end]

    f.close()
    return boundary_dict

def count_sample_variants(cm_boundaries, seg_q, vcf_file):
    sample_snp_dict = dict()

    query_boundaries = cm_boundaries[seg_q]
    start = query_boundaries[0]
    end = query_boundaries[1]
    segment_snps = []
    f = open(vcf_file, 'r')
    for line in f:
        # header
        if '#' in line:
            continue
        L = line.strip().split()
        pos = int(L[1])
        if (pos >= start and pos <= end):
            genotypes = L[9:]
            for sample in range(len(genotypes)):
                sample_gt = genotypes[sample]
                if '1' in sample_gt or '2' in sample_gt:
                    try:
                        sample_snp_dict[sample] += 1
                    except KeyError:
                        sample_snp_dict[sample] = 1
            # num_samples = count_samples(genotypes)
            # segment_snps.append(num_samples)
    f.close()
    return sample_snp_dict

def plot_hist(segment_snps, seg, chrm, cm, out_png):
    png_name = out_png + 'chr' + chrm \
               + '.seg.' + str(seg) \
               + '.variants_distribution.png'
    png_title = 'Distribution of variants for segment ' \
                + str(seg) \
                + '\n(Chromosome ' + str(chrm) + ')' \
                + '\n Segment Distance = ' + str(cm) + 'cM'

    hist_data = []
    for sample in segment_snps:
            hist_data.append(segment_snps[sample])

    plt.figure(figsize=(15, 12))
    ax1 = plt.subplot(111)
    ax1.set_title(png_title, fontsize=30)
    ax1.hist(hist_data, color='olivedrab')
    ax1.set_xlabel('number of SNPs', fontsize=20)
    # ax1.hist(segment_snps, bins=1000, color='olivedrab')
    # ax1.set_xlabel('number of samples with SNPs\n (heterozygous or homozygous alt)', fontsize=20)
    ax1.spines.right.set_visible(False)
    ax1.spines.top.set_visible(False)


    plt.savefig(png_name)

if __name__ == '__main__':
    main()

