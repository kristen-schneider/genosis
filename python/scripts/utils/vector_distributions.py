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

    vector_sizes = vector_size(segment_dir, ext)
    plot_vector_sizes(vector_sizes)

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
    num_variants = len(line)
    #for line in f:
    #    if len(line) != vector_size:
    #        print('inconsistent vector size')
    return vector_size

def get_vector_size(file):
    f = open(file, 'r')
    line = f.readline()
    vector_size = len(line)
    return vector_size

def plot_vector_sizes(vector_sizes):
    plt.figure(figsize=(15, 10))
    ax1 = plt.subplot(111)
    ax1.set_title('Segment Vector Size', fontsize=30)
    ax1.set_xlabel('vector size (number SNPs)', fontsize=20)
    ax1.hist(vector_sizes, bins=40, color='olivedrab')
    #ax1.set_xlim(xrange)
    # ax1.set_ylim(yrange)
    #ax1.set_xticks([])
    ax1.spines.right.set_visible(False)
    ax1.spines.top.set_visible(False)
    plt.savefig('vector_size.png')

if __name__ == '__main__':
    main()
