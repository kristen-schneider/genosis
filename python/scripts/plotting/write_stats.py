import argparse
import matplotlib.pyplot as plt
import statistics
import sys, os

import os, sys, inspect

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir)
sys.path.insert(0, parentdir+'/utils/')

import compute_conditions
import basic_datastructures


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--data_dir')
    parser.add_argument('--plink_ext')
    parser.add_argument('--faiss_enc_ext')
    parser.add_argument('--faiss_emb_ext')
    parser.add_argument('--faiss_idx')
    parser.add_argument('--base_name')
    parser.add_argument('--train')
    parser.add_argument('--test')
    parser.add_argument('--k')
    parser.add_argument('--num_seg')
    parser.add_argument('--out_dir')
    return parser.parse_args()

def main():
    args = get_args()
    query_IDs = basic_datastructures.get_db_q_IDs(args.test)
    database_IDs = basic_datastructures.get_db_q_IDs(args.train)

    # write average true positives and precision for each segment (encodings and embeddings)
    segment_averages = compute_segment_statistic(args, query_IDs, database_IDs, 'avg')
    write_plotting_data(args, segment_averages, 'avg', args.faiss_idx)

    # write median true positives and precision for each segment (encodings and embeddings)
    segment_medians = compute_segment_statistic(args, query_IDs, database_IDs, 'med')
    write_plotting_data(args, segment_medians, 'med', args.faiss_idx)

    return 0

def compute_segment_statistic(args, query_IDs, database_IDs, stat):
    encoding_true_positives_summary = dict()
    embedding_true_positives_summary = dict()
    encoding_precisions_summary = dict()
    embedding_precisions_summary = dict()

    for s in range(int(args.num_seg)):
        print('computing ' + stat + ' for segment ...' + str(s))
        base = args.data_dir + args.base_name
        plink_f = base + str(s) + args.plink_ext
        faiss_enc_f = base + str(s) + args.faiss_enc_ext
        faiss_emb_f = base + str(s) + args.faiss_emb_ext

        # get distances for top k
        # dicts = [queryID: [(sampleID, dist), (sampleID, dist)...], ...]

        # PLINK
        plink_dict = basic_datastructures.get_plink_distances(
            database_IDs, query_IDs, plink_f)
        plink_top = compute_conditions.compute_top_k(plink_dict, int(args.k))
        # FAISS ENCODING
        encoding_faiss_dict = basic_datastructures.get_faiss_distances(
            database_IDs, query_IDs, faiss_enc_f)
        encoding_faiss_top = compute_conditions.compute_top_k(encoding_faiss_dict, int(args.k))
        # FAISS EMBEDDING
        embedding_faiss_dict = basic_datastructures.get_faiss_distances(
            database_IDs, query_IDs, faiss_emb_f)
        embedding_faiss_top = compute_conditions.compute_top_k(embedding_faiss_dict, int(args.k))

        # compute metrics
        # get all true positives
        encoding_true_positives = compute_conditions.compute_true_positives(plink_top, encoding_faiss_top)
        embedding_true_positives = compute_conditions.compute_true_positives(plink_top, embedding_faiss_top)
        encoding_seg_TP_values = [encoding_true_positives[p] for p in query_IDs]
        embedding_seg_TP_values = [embedding_true_positives[p] for p in query_IDs]
        # get all precision
        encoding_precisions = compute_conditions.confusion_matrix(encoding_true_positives, int(args.k))
        embedding_precisions = compute_conditions.confusion_matrix(embedding_true_positives, int(args.k))
        encoding_seg_prec_values = [encoding_precisions[p] for p in query_IDs]
        embedding_seg_prec_values = [embedding_precisions[p] for p in query_IDs]
        # AVERAGE
        if stat == 'avg':
            # get avg of precision
            stat_encoding_seg_precision = (sum(encoding_seg_prec_values) / len(encoding_seg_prec_values))
            stat_embedding_seg_precision = (sum(embedding_seg_prec_values) / len(embedding_seg_prec_values))
            # get avg of true positives
            stat_encoding_seg_TP = (sum(encoding_seg_TP_values) / len(encoding_seg_TP_values))
            stat_embedding_seg_TP = (sum(embedding_seg_TP_values) / len(embedding_seg_TP_values))

        elif stat == 'med':
            # get median of precision
            stat_encoding_seg_precision = statistics.median(encoding_seg_prec_values)
            stat_embedding_seg_precision = statistics.median(embedding_seg_prec_values)
            # get median of true positives
            stat_encoding_seg_TP = statistics.median(encoding_seg_TP_values)
            stat_embedding_seg_TP = statistics.median(embedding_seg_TP_values)

        # add medians to dictionary
        encoding_true_positives_summary[s] = stat_encoding_seg_TP
        embedding_true_positives_summary[s] = stat_embedding_seg_TP
        encoding_precisions_summary[s] = stat_encoding_seg_precision
        embedding_precisions_summary[s] = stat_embedding_seg_precision

        # debugging
        #print('enc:TP:', stat_encoding_seg_TP, 'emb:TP:', stat_embedding_seg_TP,
        #      'enc:P:', stat_encoding_seg_precision, 'emb:P:', stat_embedding_seg_precision)

    stats = [encoding_true_positives_summary, embedding_true_positives_summary,
               encoding_precisions_summary, embedding_precisions_summary]
    return stats

def write_plotting_data(args, segment_stats, stat, idx):
    [ENC_TP, EMB_TP, ENC_P, EMB_P] = segment_stats
    # write tp data to outfile
    tp_enc_string = args.out_dir + stat + '.tp.encoding.' + idx
    tp_enc_data = open(tp_enc_string, 'w')
    for seg in ENC_TP.keys():
        tp_line = str(seg) + '\t' + str(ENC_TP[seg]) + '\n'
        tp_enc_data.write(tp_line)
    tp_enc_data.close()

    tp_emb_string = args.out_dir + stat + '.tp.embedding.' + idx
    tp_emb_data = open(tp_emb_string, 'w')
    for seg in EMB_TP.keys():
        tp_line = str(seg) + '\t' + str(EMB_TP[seg]) + '\n'
        tp_emb_data.write(tp_line)
    tp_emb_data.close()

    # write p data to outfile
    p_enc_string = args.out_dir + stat + '.p.encoding.' + idx
    p_enc_data = open(p_enc_string, 'w')
    for seg in ENC_P.keys():
        p_line = str(seg) + '\t' + str(ENC_P[seg]) + '\n'
        p_enc_data.write(p_line)
    p_enc_data.close()

    p_emb_string = args.out_dir + stat + '.p.embedding.' + idx
    p_emb_data = open(p_emb_string, 'w')
    for seg in EMB_P.keys():
        p_line = str(seg) + '\t' + str(EMB_P[seg]) + '\n'
        p_emb_data.write(p_line)
    p_emb_data.close()


if __name__ == '__main__':
    main()
