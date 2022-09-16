import argparse
import matplotlib.pyplot as plt


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--hd')
    parser.add_argument('--plink')
    parser.add_argument('--png')
    return parser.parse_args()


def main():
    args = get_args()
    hd_dict = read_distance_file(args.hd)
    plink_dict = read_distance_file(args.plink)
    plot_hd_plink_best_matches(hd_dict, plink_dict, args.png)


def plot_hd_plink_best_matches(hd_dict, plink_dict, png):
    # format data
    # 1. get both haplotype from ed
    hd_haplotypes = dict()
    for sample_h in hd_dict:
        plink_hap_name = sample_h.strip().split('_')[0]
        try:
            hd_haplotypes[plink_hap_name].append(hd_dict[sample_h])
        except KeyError:
            hd_haplotypes[plink_hap_name] = [hd_dict[sample_h]]

    hd_data = []
    plink_data = []
    for sample in hd_haplotypes:
        best_hd_hap = min(hd_haplotypes[sample])
        if best_hd_hap != 0.0:
            hd_data.append(best_hd_hap)
        try:
            plink_data.append(plink_dict[sample])
        except KeyError:
            print('self')
            continue

    # plot data
    plt.figure(figsize=(10, 10))
    ax1 = plt.subplot(111)
    png_title = 'Plink distance vs Hamming distance\n' \
                'for binary haplotype encodings\n' \
                'on 1 cm segments, for 1 query\n (Segment 0, CHR 8)'
    ax1.set_title(png_title)
    ax1.scatter(plink_data, hd_data, color='olivedrab')
    # format plot
    png_name = png + 'hamming_plink.png'
    ax1.set_xlabel('Plink Distance')
    ax1.set_ylabel('Hamming Distance')
    ax1.spines.top.set_visible(False)
    ax1.spines.right.set_visible(False)
    plt.savefig(png_name)

def read_distance_file(distance_file):
    sample_ID_dist_dict = dict()

    f = open(distance_file, 'r')
    query = f.readline().strip().split()[1]
    header = f.readline()
    for line in f:
        L = line.strip().split()
        sample_ID = L[0]
        dist = float(L[1])
        sample_ID_dist_dict[sample_ID] = dist
    f.close()

    return sample_ID_dist_dict



if __name__ == '__main__':
    main()

