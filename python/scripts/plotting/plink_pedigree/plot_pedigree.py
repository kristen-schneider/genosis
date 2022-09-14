import argparse
from collections import defaultdict
import matplotlib.pyplot as plt
import numpy as np

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--plink')
    parser.add_argument('--samples')
    parser.add_argument('--png')
    return parser.parse_args()

def main():
    args = get_args()
    samples = read_samples(args.samples)
    plink_dict = read_plink_dict(args.plink)
    plot_data = make_plot_data(samples, plink_dict)
    plot_plink_histogram(plot_data, samples, args.png)

def plot_plink_histogram(plot_data_dict, samples, png):
    plot_data_list = [plot_data_dict[sample] for sample in plot_data_dict]

    png_name = png + 'relatedness.png'
    #plt.figure(figsize=(28, 15))
    ax1 = plt.subplot(111)
    # im = sns.heatmap(plot_data_list, cmap='mako')
    png_title = 'sample relatedness'
    ax1.set_title(png_title)
    im = ax1.imshow(plot_data_list, cmap='inferno', vmin=0.95)
    cbar = plt.colorbar(im, extend='min')
    ax1.set_xticks(np.arange(len(samples)), labels=samples, rotation=45)
    ax1.set_yticks(np.arange(len(samples)), labels=samples) 
    plt.savefig(png_name)

def make_plot_data(samples, plink_dict):
    plot_dict = dict()
    for sample_1 in samples:
        curr_sample = plink_dict[sample_1]
        try:
            plot_dict[sample_1] = [curr_sample[sample_2] for sample_2 in samples]
        except KeyError:
            print(sample_1)
    return plot_dict


def read_plink_dict(plink_file):
    plink_dict = defaultdict(dict)
    f = open(plink_file, 'r')
    for h in range(3):
        f.readline()
    for line in f:
        L = line.strip().split()
        SID1 = L[0]
        SID2 = L[1]
        dist = float(L[2])
        # add self to dictionary
        plink_dict[SID1][SID1] = 0.0
        plink_dict[SID1][SID2] = dist
        plink_dict[SID2][SID1] = dist

    f.close()
    return plink_dict

def read_samples(samples_file):
    samples_list = []
    f = open(samples_file, 'r')
    for line in f:
        L = line.strip().split()
        SID = L[0]
        samples_list.append(SID)
    f.close()
    return samples_list


if __name__ == '__main__':
    main()
