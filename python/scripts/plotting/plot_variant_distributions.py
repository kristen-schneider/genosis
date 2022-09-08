import argparse
import matplotlib.pyplot as plt
import sys, os

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--vcf')
    parser.add_argument('--data_dir')
    parser.add_argument('--out_png')
    parser.add_argument('--cm')
    parser.add_argument('--chr')
    parser.add_argument('--seg')
    return parser.parse_args()

def main():
    args = get_args()
    boundary_file = args.data_dir \
                    + 'chr' + str(args.chr) \
                    + '.' + str(args.cm) \
                    + 'cm-segment-boundaries'
    max_snps = dict()
    cm_boundaries = read_boundary_data(boundary_file)
    num_segs = len(cm_boundaries)
    for s in range(0, num_segs, 50):
        print('segment ', s)
        sample_snps_dict = count_sample_variants(cm_boundaries, s, args.vcf)
        segment_counts = list(sample_snps_dict.values())
        max_snps[s] = max(segment_counts)
        #plot_hist(sample_snps_dict, s, args.chr, args.cm, args.out_png)
        # print('counting samples for segment ', s)
        # segment_snps = read_vcf(cm_boundaries, s, args.vcf)
        # plot_hist(segment_snps, s, args.chr, args.cm, args.out_png)
    
    for seg in max_snps:
        print (seg, ', ', max_snps[seg])

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



def read_boundary_data(boundary_file):
    boundary_dict = dict()
    f = open(boundary_file)

    header = None
    for line in f:
        if header == None:
            header = line
            continue
        L = line.strip().split(',')
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
    png_title = 'Frequency of variants for segment ' \
                + str(seg) \
                + '\n(Chromosome ' + str(chrm) + ')' \
                + '\n Segment Distance = ' + str(cm) + 'cM'

    hist_data = []
    for sample in segment_snps:
        hist_data.append(segment_snps[sample])

    plt.figure(figsize=(15, 12))
    ax1 = plt.subplot(111)
    ax1.set_title(png_title, fontsize=30)
    ax1.hist(hist_data, bins=30, color='olivedrab')
    ax1.set_xlabel('number of SNPs', fontsize=20)
    ax1.set_ylabel('number of samples', fontsize=20)
    # ax1.hist(segment_snps, bins=1000, color='olivedrab')
    # ax1.set_xlabel('number of samples with SNPs\n (heterozygous or homozygous alt)', fontsize=20)
    ax1.spines.right.set_visible(False)
    ax1.spines.top.set_visible(False)

    plt.savefig(png_name)

if __name__ == '__main__':
    main()

