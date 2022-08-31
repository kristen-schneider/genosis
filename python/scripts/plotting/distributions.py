import argparse
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt

import os, sys, inspect

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir)
sys.path.insert(0, parentdir+'/utils/')

import basic_datastructures

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--data_dir')
    parser.add_argument('--plink_ext')
    parser.add_argument('--faiss_enc_ext')
    parser.add_argument('--faiss_emb_ext')
    parser.add_argument('--train')
    parser.add_argument('--test')
    parser.add_argument('--num_seg')
    parser.add_argument('--out_dir')
    return parser.parse_args()

def main():
    args = get_args()
    query_IDs = basic_datastructures.get_db_q_IDs(args.test)
    database_IDs = basic_datastructures.get_db_q_IDs(args.train)

    plink_distributions(database_IDs, query_IDs, args.data_dir, int(args.num_seg), args.plink_ext, args.out_dir + 'plink_distribution.png')
    # faiss_distributions(database_IDs, query_IDs, args.data_dir, int(args.num_seg), args.faiss_enc_ext, args.out_dir + 'faiss_enc_distribution.png')
    # faiss_distributions(database_IDs, query_IDs, args.data_dir, int(args.num_seg), args.faiss_emb_ext, args.out_dir + 'faiss_emb_distribution.png')

def plink_distributions(database_IDs, query_IDs, data_dir, num_seg,
                        plink_ext, out_png):
    """
    plots boxplot of distances
        reported by plink for many segments
    """

    plink_seg_data = dict()
    for s in range(int(num_seg)):
        print('getting distances for segment...' + str(s))

        for segment_file in os.listdir(data_dir):
            f = os.path.join(data_dir, segment_file)
            base = 'seg.' + str(s)
            if (base in f and plink_ext in f):
                plink_dict = basic_datastructures.get_plink_distances(database_IDs, query_IDs, f)
                plink_data = []
                for q in query_IDs:
                    for pair in plink_dict[q]:
                        plink_data.append(pair[1])

                plink_seg_data[s] = plink_data

    faiss_plot_data = []
    for seg in range(num_seg):
        try:
            faiss_plot_data.append(plink_seg_data[seg])
        except KeyError:
            continue

    # plotting
    plt.figure(figsize=(15, 15))
    # plot for many queries in one segment file

    ax = plt.subplot(111)
    ax.boxplot(faiss_plot_data, showfliers=False)
    ax.set_xticks(range(num_seg), rotation=90)
    ax.spines.right.set_visible(False)
    ax.spines.top.set_visible(False)

    plt.title('Distribution of Plink Distances')
    plt.xlabel('Segment')
    plt.ylabel('Plink Distance distance')
    plt.savefig(out_png)

def faiss_distributions(database_IDs, query_IDs, data_dir, num_seg,
                        faiss_ext, out_png):
    """
    plots boxplot of distances
    reported by faiss for many segments
    :param faiss_file:
    :return:
    """

    faiss_seg_data = dict()
    for s in range(int(num_seg)):
        print('getting distances for segment...' + str(s))

        for segment_file in os.listdir(data_dir):
            f = os.path.join(data_dir, segment_file)
            base = 'seg.' + str(s)
            if (base in f and faiss_ext in f):
                faiss_dict = basic_datastructures.get_faiss_distances(database_IDs, query_IDs, f)
                faiss_data = []
                for q in query_IDs:
                    for pair in faiss_dict[q]:
                        if pair[1] > 5:
                            print('segment: ', s, 'query: ', q, pair)
                        faiss_data.append(pair[1])

                faiss_seg_data[s] = faiss_data

    faiss_plot_data = []
    for seg in range(num_seg):
        try: faiss_plot_data.append(faiss_seg_data[seg])
        except KeyError: continue

    # plotting
    plt.figure(figsize=(15, 15))
    # plot for many queries in one segment file

    ax = plt.subplot(111)
    ax.boxplot(faiss_plot_data, showfliers=False)
    ax.set_xticks(range(num_seg), rotation=90)
    ax.spines.right.set_visible(False)
    ax.spines.top.set_visible(False)

    plt.title('Distribution of FAISS Distances')
    plt.xlabel('Segment')
    plt.ylabel('FAISS distance')
    plt.savefig(out_png)

if __name__ == '__main__':
    main()
