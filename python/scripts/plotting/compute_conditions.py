import sys
import os, sys, inspect

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir)

from utils import basic_datastructures

def compute_false_positives(plink_top, faiss_top):
    query_false_positives = dict()

    for query_ID in plink_top.keys():
        query_FP_count = 0
        query_top_IDs = []

        plink_query = plink_top[query_ID]
        faiss_query = faiss_top[query_ID]

        for plink_top_hit in plink_query:
            query_top_IDs.append(plink_top_hit[0])
        for faiss_top_hit in faiss_query:
            faiss_ID = faiss_top_hit[0]
            if faiss_ID not in query_top_IDs:
                query_FP_count += 1

        query_false_positives[query_ID] = query_FP_count
    return query_false_positives

def compute_true_positives(plink_top, faiss_top):
    query_true_positives = dict()

    for query_ID in plink_top.keys():
        query_TP_count = 0
        query_top_IDs = []

        plink_query = plink_top[query_ID]
        faiss_query = faiss_top[query_ID]

        for plink_top_hit in plink_query:
            query_top_IDs.append(plink_top_hit[0])
        for faiss_top_hit in faiss_query:
            faiss_ID = faiss_top_hit[0]
            if faiss_ID in query_top_IDs:
                query_TP_count += 1

        query_true_positives[query_ID] = query_TP_count
    return query_true_positives


def compute_positives(plink_top, faiss_top):
    '''
    Returns the top k hits for all queries
    :param plink_dict: sorted dictionary {queryID: [(databaseID, dist), (databaseID, dist), ...]}
    :param faiss_dict: file with IDs of all
    :param k: top k
    :return: number of unique positives found between
    both methods (plink, faiss) for each query
    '''
    query_positives = dict()

    for query_ID in plink_top.keys():
        query_top_IDs = []

        plink_query = plink_top[query_ID]
        faiss_query = faiss_top[query_ID]

        for plink_top_hit in plink_query:
            query_top_IDs.append(plink_top_hit[0])
        for faiss_top_hit in faiss_query:
            faiss_ID = faiss_top_hit[0]
            if faiss_ID not in query_top_IDs:
                query_top_IDs.append(faiss_ID)

        query_positives[query_ID] = len(query_top_IDs)
    return query_positives

def compute_top_k(data_dict, k):
    all_query_top = dict()
    for q in range(len(data_dict)):
        query_top = []
        dict_keys = list(data_dict.keys())
        data_q = dict_keys[q]
        topK_hits = data_dict[data_q][:k]
        for tuple in topK_hits:
            query_top.append(tuple)
        all_query_top[data_q] = query_top
    return all_query_top




