import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt

import sys

from python.utils import basic_datastructures
# sys.path.insert(1, '/home/sdp/genotype-encoding/python/utils/')
# import basic_datastructures


def main():
    # system arguments
    sample_ID_file = sys.argv[1]
    plink_file = sys.argv[2]
    faiss_encoding_file = sys.argv[3]
    faiss_embedding_file = sys.argv[4]
    plink_plot = sys.argv[5]
    faiss_encoding_plot = sys.argv[6]
    faiss_embedding_plot = sys.argv[7]
    num_queries = 10

    sample_ID_list = basic_datastructures.get_sample_ID_list(sample_ID_file)

    # plink_distributions(plink_file, plink_plot)
    # faiss_distributions(sample_ID_list, faiss_encoding_file, faiss_encoding_plot, num_queries)
    faiss_distributions(sample_ID_list, faiss_embedding_file, faiss_embedding_plot, num_queries)

def plink_distributions(plink_file, plink_plot):
    """
    plots histogram of distances
    reported by plink for one plink file
    :param plink_file:
    :return:
    """
    plink_distances = []
    f = open(plink_file, 'r')
    header = f.readline()
    for line in f:
        l = line.strip().split()
        dist = float(l[11])
        plink_distances.append(dist)
    f.close()

    # plotting
    plt.figure(figsize=(15, 15))
    ax = plt.subplot(111)
    ax.hist(plink_distances, bins=30, color='salmon')
    ax.spines.right.set_visible(False)
    ax.spines.top.set_visible(False)
    plt.title('Distribution of Plink Distances')
    plt.xlabel('Plink Distance')
    plt.ylabel('Frequency')
    plt.xlim(0.944, 1)
    # plt.xticks(rotation=90)
    plt.savefig(plink_plot)

def faiss_distributions(sample_ID_list, faiss_file, faiss_plot, num_queries):
    """
    plots histogram of distances
    reported by faiss for one faiss file
    :param faiss_file:
    :return:
    """
    faiss_dict = basic_datastructures.get_faiss_distances_dict(sample_ID_list, faiss_file)
    query_IDs = sample_ID_list[:num_queries]
    faiss_distances = []

    # plotting
    plt.figure(figsize=(15, 15))

    colors = cm.rainbow(np.linspace(0, 1, num_queries))
    # plot for many queries in one segment file
    for q, c in zip(query_IDs, colors):
        query_ID = q
        print(query_ID)
        for id in sample_ID_list:
            curr_distance = faiss_dict[query_ID][id]
            faiss_distances.append(curr_distance)
        ax = plt.subplot(111)
        ax.hist(faiss_distances, bins=30, color=c, alpha=0.2)
        ax.spines.right.set_visible(False)
        ax.spines.top.set_visible(False)

    plt.title('Distribution of FAISS Distances')
    plt.xlabel('FAISS Distance')
    plt.ylabel('Frequency')
    plt.savefig(faiss_plot)

if __name__ == '__main__':
    main()
