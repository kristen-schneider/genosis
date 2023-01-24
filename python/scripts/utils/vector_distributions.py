import argparse
import matplotlib.pyplot as plt
import os

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--segment_dir')
    parser.add_argument('--ext')
    return parser.parse_args()

def main():
    args = get_args()
    segment_dir = args.segment_dir
    ext = args.ext

    #vector_sizes = vector_size(segment_dir, ext)
    #plot_vector_sizes(vector_sizes)

    max_variants = num_variants(segment_dir, ext)
    plot_vector_sizes(max_variants)

def num_variants(segment_dir, ext):
    max_variants = []

    for file in os.listdir(segment_dir):
        if ext in file:
            print(file)
            file_max_size = get_num_variants(segment_dir+file)
            max_variants.append(file_max_size)
    return max_variants

def vector_size(segment_dir, ext):
    vector_sizes = []

    for file in os.listdir(segment_dir):
        if ext in file:
            #print(file)
            file_vec_size = get_vector_size(segment_dir+file)
            vector_sizes.append(file_vec_size)
    return vector_sizes

def get_num_variants(file):
    f = open(file, 'r')
    line = f.readline()
    max_variants = -1
    for line in f:
        num_variants = 0
        for v in line:
            if v == '1':
                num_variants += 1
        if num_variants > max_variants:
            max_variants = num_variants
    return max_variants

def get_vector_size(file):
    f = open(file, 'r')
    line = f.readline()
    vector_size = len(line)
    return vector_size

def plot_vector_sizes(vector_sizes):
    plt.figure(figsize=(15, 10))
    ax1 = plt.subplot(111)
    ax1.set_title('Segment Number Variants', fontsize=30)
    ax1.set_xlabel('Number Variants (number SNPs)', fontsize=20)
    ax1.hist(vector_sizes, bins=40, color='steelblue')
    #ax1.set_xlim(xrange)
    # ax1.set_ylim(yrange)
    #ax1.set_xticks([])
    ax1.spines.right.set_visible(False)
    ax1.spines.top.set_visible(False)
    plt.savefig('num_variants.png')

if __name__ == '__main__':
    main()
