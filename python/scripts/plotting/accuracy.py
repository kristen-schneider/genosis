import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt

import os, sys, inspect

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir)

import basic_datastructures

def main():
    # system arguments
    sample_ID_file = sys.argv[1]
    plink_file = sys.argv[2]
    faiss_encoding_file = sys.argv[3]
    faiss_embedding_file = sys.argv[4]
    faiss_ranks_enc_plot = sys.argv[5]
    faiss_ranks_emb_plot = sys.argv[6]


    plink_upper_cutoff = 0.985
    plink_lower_cutoff = 0.955
    sample_ID_list = basic_datastructures.get_sample_ID_list(sample_ID_file)

    positive_sample_pairs, negative_sample_pairs = \
        get_pos_neg_sample_pairs(plink_file, plink_upper_cutoff, plink_lower_cutoff)

    # encodings
    faiss_ranks_dict_enc = basic_datastructures.get_faiss_rankings_dict(sample_ID_list, faiss_encoding_file)
    print('positive encodings')
    faiss_pos_ranks_enc = get_ranks_from_pairs(positive_sample_pairs, faiss_ranks_dict_enc)
    print('negative encodings')
    faiss_neg_ranks_enc = get_ranks_from_pairs(negative_sample_pairs, faiss_ranks_dict_enc)
    plot_rank_frequency(faiss_pos_ranks_enc, faiss_neg_ranks_enc, faiss_ranks_enc_plot)

    # embeddings
    faiss_ranks_dict_emb = basic_datastructures.get_faiss_rankings_dict(sample_ID_list, faiss_embedding_file)
    print('positive embeddings')
    faiss_pos_ranks_emb = get_ranks_from_pairs(positive_sample_pairs, faiss_ranks_dict_emb)
    print('negative embeddings')
    faiss_neg_ranks_emb = get_ranks_from_pairs(negative_sample_pairs, faiss_ranks_dict_emb)
    plot_rank_frequency(faiss_pos_ranks_emb, faiss_neg_ranks_emb, faiss_ranks_emb_plot)


    debug=''

def plot_rank_frequency(faiss_pos_ranks, faiss_neg_ranks, faiss_ranks_plot):

    # plotting
    plt.figure(figsize=(15, 15))
    ax = plt.subplot(111)
    ax.hist(faiss_pos_ranks, bins=100, color='yellowgreen', alpha=0.7)
    ax.hist(faiss_neg_ranks, bins=100, color='salmon', alpha=0.7)
    ax.spines.right.set_visible(False)
    ax.spines.top.set_visible(False)
    plt.title('Distribution of FAISS Ranks for Positive Pairs')
    plt.xlabel('FAISS Rank')
    plt.ylabel('Frequency')
    plt.xlim(0, 2548)
    # plt.xticks(rotation=90)
    plt.savefig(faiss_ranks_plot)

def get_pos_neg_sample_pairs(plink_file, plink_upper_cutoff, plink_lower_cutoff):
    positive_sample_pairs = []
    negative_sample_pairs = []

    f = open(plink_file, 'r')
    header = f.readline()
    for line in f:
        l = line.strip().split()
        sample1 = l[1]
        sample2 = l[3]
        distance = float(l[11])
        pair = (sample1, sample2)
        if distance >= plink_upper_cutoff:
            positive_sample_pairs.append(pair)
        elif distance <= plink_lower_cutoff:
            negative_sample_pairs.append(pair)

    f.close()
    return positive_sample_pairs, negative_sample_pairs

def get_ranks_from_pairs(pairs, faiss_ranks_dict):
    """
    given a list of pairs which meet the
    distance cutoff requirement in plink,
    find their corresponding faiss rank

    :param positive_pairs:
    :param faiss_ranks_dict:
    :return:
    """
    faiss_ranks = []
    top_ten = 0

    for pair in pairs:
        sample1 = pair[0]
        sample2 = pair[1]

        pair_faiss_rank = faiss_ranks_dict[sample1][sample2]
        faiss_ranks.append(pair_faiss_rank)
    return faiss_ranks

if __name__ == '__main__':
    main()
