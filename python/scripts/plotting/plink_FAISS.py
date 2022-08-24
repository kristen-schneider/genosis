import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from scipy.stats import kendalltau
from scipy.stats import mannwhitneyu

import os, sys, inspect

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir)

import basic_datastructures

sample_ID_file = sys.argv[1]
plink_file = sys.argv[2]
FAISS_file = sys.argv[3]
plink_euclidean_distances_png = sys.argv[4]
plink_euclidean_rank_png = sys.argv[5]


def main():
    # x = np.random.rand(100)
    # y = np.random.rand(100)
    # U1, p = mannwhitneyu(x, y, method="asymptotic")
    # print(U1, p)
    # U2 = (len(x) * len(x)) - U1
    # print(U2)

    sampleIDs = basic_datastructures.get_sample_ID_list(sample_ID_file)
    query_samples = range(1) # number of queries to be plotted
    topX = range(10, 300) # number of top queries to plot

    # get plink data structures
    plink_dict = basic_datastructures.get_plink_dict(plink_file)

    # get FAISS data structures
    faiss_distances = basic_datastructures.get_faiss_distances_dict(sampleIDs, FAISS_file)
    faiss_ranks = basic_datastructures.get_faiss_rankings(sampleIDs, FAISS_file)

    plot_distances(sampleIDs, query_samples, plink_dict, faiss_distances)
    # plot_recall_percent(sampleIDs, query_samples, plink_dict, faiss_ranks, topX)
    # plot_rank_compare_metric(sampleIDs, query_samples, plink_dict, faiss_ranks, topX)

    x = 'debug'

def plot_distances(sampleIDs, query_samples, plink_dict, faiss_distances):

    # plot
    plt.figure(figsize=(15, 15))
    colors = cm.rainbow(np.linspace(0, 1, len(query_samples)))
    # plot for many queries in one segment file
    for q, c in zip(query_samples, colors):
        query_ID = sampleIDs[q]
        print(query_ID)

        [plink_x, faiss_y] = get_query_XY_distances(query_ID, sampleIDs, plink_dict, faiss_distances)
        ax = plt.subplot(111)
        ax.scatter(plink_x, faiss_y, color=c)
        ax.spines.right.set_visible(False)
        ax.spines.top.set_visible(False)

    # faiss distances on the x axis
    # plink distances on y axis
    plt.title('FAISS L2 distances vs PLINK distances')
    plt.xlabel('PLINK IBD distance')
    plt.ylabel('FAISS L2 similarity')
    plt.scatter(plink_x, faiss_y, color='red')
    plt.savefig(plink_euclidean_distances_png)

def plot_recall_percent(sampleIDs, query_samples, plink_dict, faiss_ranks, topX):
    # plot
    # plt.figure(figsize=(15, 15))
    colors = cm.rainbow(np.linspace(0, 1, len(query_samples)))
    query_ID = sampleIDs[0] # one segment, one query,

    recall_percent_list = []
    for x in topX:
        [plink_x, faiss_y] = get_query_XY_ranks(query_ID, sampleIDs,
                                                plink_dict, faiss_ranks)
        [plink_x, faiss_y] = [plink_x[:x], faiss_y[:x]]
        recall_percent = get_percent_recall(plink_x, faiss_y)
        recall_percent_list.append(recall_percent)

    ax = plt.subplot(111)
    ax.scatter(topX, recall_percent_list)
    ax.spines.right.set_visible(False)
    ax.spines.top.set_visible(False)

    plt.title('FAISS L2 recalled for top X')
    plt.xlabel('X')
    plt.ylabel('Recall Percent')
    plt.ylim([0, 1.1])
    plt.savefig(plink_euclidean_rank_png)

def plot_rank_compare_metric(sampleIDs, query_samples, plink_dict, faiss_ranks, topX):
    # plot
    # plt.figure(figsize=(15, 15))
    query_ID = sampleIDs[0] # one segment, one query,

    rank_compare_correlation_list = []
    for x in topX:
        [plink_x, faiss_y] = get_query_XY_ranks(query_ID, sampleIDs,
                                                plink_dict, faiss_ranks)
        [plink_x, faiss_y] = [plink_x[:x], faiss_y[:x]]
        # tau, pval = kendalltau(np.array(plink_x), np.array(faiss_y))
        U1, pval = mannwhitneyu(np.array(plink_x), np.array(faiss_y),
                                method='asymptotic')

        rank_compare_correlation_list.append(pval)

    ax = plt.subplot(111)
    ax.scatter(topX, rank_compare_correlation_list)
    ax.spines.right.set_visible(False)
    ax.spines.top.set_visible(False)

    plt.title('FAISS L2 MannWhitneyU P-Value for top X')
    plt.xlabel('X')
    plt.ylabel('MannWhitneyU P-Value')
    # plt.ylim([0, 1.1])
    plt.savefig(plink_euclidean_rank_png)

def get_query_XY_distances(query_ID, sampleIDs, plink_dict, faiss_distances):
    plink_query = plink_dict[query_ID]
    plink_x = []
    for sample2 in sampleIDs:
        try:
            plink_x.append(plink_query[sample2])
        except KeyError:
            continue

    faiss_query = faiss_distances[query_ID]
    faiss_y = []
    for sample2 in sampleIDs:
        # faiss will match self to self
        if(sample2 != query_ID):
            try:
                faiss_y.append(faiss_query[sample2])
            except KeyError:
                continue

    return[plink_x, faiss_y]

def get_query_XY_ranks(query_ID, sampleIDs, plink_dict, faiss_ranks):
    plink_query = plink_dict[query_ID]
    plink_temp = dict(sorted(plink_query.items(), key=lambda item: item[1], reverse=True))
    plink_x = []
    for k in list(plink_temp.keys()):
        plink_x.append(sampleIDs.index(k))


    faiss_query = faiss_ranks[query_ID]
    faiss_y = []
    for rank in faiss_query:
        sample2 = sampleIDs[rank]
        # faiss will match self to self,
        # we need to ignore this because plink does not
        if(sample2 != query_ID):
            try:
                faiss_y.append(rank)
            except KeyError:
                continue

    return[plink_x, faiss_y]

def get_percent_recall(plink_x, faiss_y):
    count_match = 0
    for x in plink_x:
        if x in faiss_y:
            count_match += 1
    return (count_match/len(plink_x))

if __name__ == '__main__':
    main()
