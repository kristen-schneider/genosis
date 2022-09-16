import argparse
import matplotlib.pyplot as plt

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--ed')
    parser.add_argument('--hd')
    parser.add_argument('--png')
    return parser.parse_args()

def main():
    args = get_args()
    ed_dict = read_distance_file(args.ed)
    hd_dict = read_distance_file(args.hd)
    plot_ed_hd(ed_dict, hd_dict, args.png)

def plot_ed_hd(ed_dict, hd_dict, png):
    # format data
    ed_data = []
    hd_data = []
    for sample_h in ed_dict:
        ed_data.append(ed_dict[sample_h])
        hd_data.append(hd_dict[sample_h])

    # plot data
    plt.figure(figsize=(10, 10))
    ax1 = plt.subplot(111)
    png_title = 'Euclidean distance vs Hamming distance\n' \
                'for binary haplotype encodings\n' \
                'on 1 cm segments for 1 query\n (Segment 0, CHR 8)'
    ax1.set_title(png_title)
    ax1.scatter(ed_data, hd_data, color='olivedrab')
    # format plot
    png_name = png + 'euclidean_hamming.png'
    ax1.set_ylabel('Hamming Distance')
    ax1.set_xlabel('Euclidean Distance')
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
