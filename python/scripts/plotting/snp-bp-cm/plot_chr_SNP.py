import argparse
import matplotlib.pyplot as plt

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--data')
    parser.add_argument('--png')
    parser.add_argument('--cm')
    parser.add_argument('--chrm')

    return parser.parse_args()

def main():
    args = get_args()
    data_dict = read_data_file(args.data)
    plot_snps(data_dict, args.png)

def plot_snps(data_dict, png_name):
    data_list = []
    for seg in data_dict:
        data_list.append(data_dict[seg])

    # plt.figure(figsize=(28, 15))
    ax1 = plt.subplot(111)
    png_title = 'Frequency of maximum SNP counts\n' \
                ' for 1 cm segments on CHR 8'
    ax1.set_title(png_title)
    ax1.hist(data_list, bins=25, color='olivedrab')
    ax1.set_xlim([0, 12000])
    ax1.set_xlabel('Max number of SNPs in a sample')
    ax1.spines.top.set_visible(False)
    ax1.spines.right.set_visible(False)
    plt.savefig(png_name)


def read_data_file(data_file):
    data_dict = dict()
    f = open(data_file, 'r')
    header = f.readline()
    for line in f:
        L = line.strip().split()
        seg = L[0]
        SNPs = int(L[1])
        data_dict[seg] = SNPs
    f.close()
    return data_dict

if __name__ == '__main__':
    main()
