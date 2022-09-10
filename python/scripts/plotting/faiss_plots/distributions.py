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
    parser.add_argument('--base_name')
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

    plink_distributions(database_IDs, query_IDs, args)
    faiss_distributions(database_IDs, query_IDs, args, 'enc')
    faiss_distributions(database_IDs, query_IDs, args, 'emb')

def plink_distributions(database_IDs, query_IDs, args):
    """
    plots boxplot of distances
        reported by plink for many segments
    """
    print('PLINK')
    plink_seg_data = dict()
    for s in range(int(args.num_seg)):
        print('getting distances for segment ...' + str(s))
        base = args.data_dir + args.base_name
        plink_f = base + str(s) + args.plink_ext
        plink_dict = basic_datastructures.get_plink_distances(
            database_IDs, query_IDs, plink_f)
        plink_data = []
        for q in query_IDs:
            for pair in plink_dict[q]:
                plink_data.append(pair[1])

        plink_seg_data[s] = plink_data

    # Get in list of lists for plotting
    plot_data = []
    for seg in range(int(args.num_seg)):
        try:
            plot_data.append(plink_seg_data[seg])
        except KeyError:
            continue

    # plotting
    plt.figure(figsize=(30, 15))
    # plot for many queries in one segment file

    ax = plt.subplot(111)
    ax.boxplot(plot_data, showfliers=False)
    ax.spines.right.set_visible(False)
    ax.spines.top.set_visible(False)

    plt.title('Distribution of Plink Distances')
    plt.xlabel('Segment')
    plt.ylabel('Plink Distance distance')
    plink_png = args.out_dir + 'plink.distributions.png'
    plt.savefig(plink_png)

def faiss_distributions(database_IDs, query_IDs, args, vector_type):
    """
    plots boxplot of distances
    reported by faiss for many segments
    :param faiss_file:
    :return:
    """
    print('FAISS')
    faiss_seg_data = dict()
    for s in range(int(args.num_seg)):
        print('getting distances for segment ...' + str(s))
        base = args.data_dir + args.base_name
        plink_f = base + str(s) + args.plink_ext
        if vector_type == 'enc':
            faiss_f = base + str(s) + args.faiss_enc_ext
        elif vector_type == 'emb':
            faiss_emb_f = base + str(s) + args.faiss_emb_ext
        # FAISS
        faiss_dict = basic_datastructures.get_faiss_distances(
            database_IDs, query_IDs, plink_f)
        faiss_data = []
        for q in query_IDs:
            try:
                for pair in faiss_dict[q]:
                    faiss_data.append(pair[1])
            except KeyError:
                continue

        faiss_seg_data[s] = faiss_data

    # Get in list of lists for plotting
    plot_data = []
    for seg in range(int(args.num_seg)):
        try: plot_data.append(faiss_seg_data[seg])
        except KeyError: continue

    # plotting
    plt.figure(figsize=(30, 15))
    # plot for many queries in one segment file

    ax = plt.subplot(111)
    ax.boxplot(plot_data, showfliers=False)
    ax.spines.right.set_visible(False)
    ax.spines.top.set_visible(False)

    plt.title('Distribution of FAISS Distances')
    plt.xlabel('Segment')
    plt.ylabel('FAISS distance')
    if vector_type == 'enc':
        faiss_png = args.out_dir + 'enc.faiss.distributions.png'
    elif vector_type == 'emb':
        faiss_png = args.out_dir + 'emb.faiss.distributions.png'
    plt.savefig(faiss_png)

if __name__ == '__main__':
    main()
