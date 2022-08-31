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
    parser.add_argument('--k')
    parser.add_argument('--train')
    parser.add_argument('--test')
    parser.add_argument('--out_dir')
    parser.add_argument('--num_seg')
    parser.add_argument('--faiss_index')
    return parser.parse_args()

def main():
    args = get_args()
    query_IDs = basic_datastructures.get_db_q_IDs(args.test)
    database_IDs = basic_datastructures.get_db_q_IDs(args.train)

    # write average true positives and precision for each segment (encodings and embeddings)
    segment_averages = compute_average_segment_statistic(args, query_IDs, database_IDs)
    write_plotting_data(args, segment_averages)

    # write median true positives and precision for each segment (encodings and embeddings)
    segment_medians = compute_median_segment_statistic(args, query_IDs, database_IDs)
    write_plotting_data(args, segment_medians)

    # write 

    return 0

def write_plotting_data(args, segment_stats, stat):
    [ENC_TP, EMB_TP, ENC_P, EMB_P] = segment_stats
    # write tp data to outfile
    tp_enc_string = args.data_dir + stat + '.tp.encoding.' + args.faiss_index
    tp_enc_data = open(tp_enc_string, 'w')
    for seg in ENC_TP.keys():
        tp_line = str(seg) + '\t' + str(ENC_TP[seg]) + '\n'
        tp_enc_data.write(tp_line)
    tp_enc_data.close()

    tp_emb_string = args.data_dir + stat + '.tp.embedding.' + args.faiss_index
    tp_emb_data = open(tp_emb_string, 'w')
    for seg in EMB_TP.keys():
        tp_line = str(seg) + '\t' + str(EMB_TP[seg]) + '\n'
        tp_emb_data.write(tp_line)
    tp_emb_data.close()

    # write p data to outfile
    p_enc_string = args.data_dir + stat + '.p.encoding.' + args.faiss_index
    p_enc_data = open(p_enc_string, 'w')
    for seg in ENC_P.keys():
        p_line = str(seg) + '\t' + str(ENC_P[seg]) + '\n'
        p_enc_data.write(p_line)
    p_enc_data.close()

    p_emb_string = args.data_dir + stat + '.p.embedding.' + args.faiss_index
    p_emb_data = open(p_emb_string, 'w')
    for seg in EMB_P.keys():
        p_line = str(seg) + '\t' + str(EMB_P[seg]) + '\n'
        p_emb_data.write(p_line)
    p_emb_data.close()

def compute_median_segment_statistic(args, query_IDs, database_IDs):
    # key = segment idx, value = median value (tp, p)
    md_encoding_true_positives = dict()
    md_embedding_true_positives = dict()
    md_encoding_precisions = dict()
    md_embedding_precisions = dict()

    for s in range(int(args.num_seg)):
        print('computing metrics for segment ...' + str(s))

        # get distances for top k
        for segment_file in os.listdir(args.data_dir):
            f = os.path.join(args.data_dir, segment_file)
            base = 'seg.' + str(s)
            # segment's plink file
            if (base in f and args.plink_ext in f):
                s_plink_file = f
                # plink dictionary
                plink_dict = basic_datastructures.get_plink_distances(
                    database_IDs, query_IDs, s_plink_file)
                plink_top = compute_conditions.compute_top_k(plink_dict, int(args.k))

            # segment's faiss file for embeddings
            elif (base in f and args.faiss_enc_ext in f):
                s_faiss_enc_file = f
                encoding_faiss_dict = basic_datastructures.get_faiss_distances(
                    database_IDs, query_IDs, s_faiss_enc_file)
                encoding_faiss_top = compute_conditions.compute_top_k(encoding_faiss_dict, int(args.k))

            # segment's faiss file for encodings
            elif (base in f and args.faiss_emb_ext in f):
                s_faiss_emb_file = f
                embedding_faiss_dict = basic_datastructures.get_faiss_distances(
                    database_IDs, query_IDs, s_faiss_emb_file)
                embedding_faiss_top = compute_conditions.compute_top_k(embedding_faiss_dict, int(args.k))

        # compute metrics
        # get all true positives
        encoding_true_positives = compute_conditions.compute_true_positives(plink_top, encoding_faiss_top)
        embedding_true_positives = compute_conditions.compute_true_positives(plink_top, embedding_faiss_top)
        encoding_seg_TP_values = [encoding_true_positives[p] for p in query_IDs]
        embedding_seg_TP_values = [embedding_true_positives[p] for p in query_IDs]
        # get avg of true positives
        avg_encoding_seg_TP = (sum(encoding_seg_TP_values) / len(encoding_seg_TP_values))
        avg_embedding_seg_TP = (sum(embedding_seg_TP_values) / len(embedding_seg_TP_values))
        # get median of true positives
        med_encoding_seg_TP = statistics.median(encoding_seg_TP_values)
        med_embedding_seg_TP = statistics.median(embedding_seg_TP_values)

        # get all precision
        encoding_precisions = compute_conditions.confusion_matrix(encoding_true_positives, int(args.k))
        embedding_precisions = compute_conditions.confusion_matrix(embedding_true_positives, int(args.k))
        encoding_seg_prec_values = [encoding_precisions[p] for p in query_IDs]
        embedding_seg_prec_values = [embedding_precisions[p] for p in query_IDs]
        # get avg of precision
        avg_encoding_seg_precision = (sum(encoding_seg_prec_values) / len(encoding_seg_prec_values))
        avg_embedding_seg_precision = (sum(embedding_seg_prec_values) / len(embedding_seg_prec_values))
        # get median of precision
        med_encoding_seg_precision = statistics.median(encoding_seg_prec_values)
        med_embedding_seg_precision = statistics.median(embedding_seg_prec_values)

        # add medians to dictionary
        md_encoding_true_positives[s] = med_encoding_seg_TP
        md_embedding_true_positives[s] = med_embedding_seg_TP
        md_encoding_precisions[s] = med_encoding_seg_precision
        md_embedding_precisions[s] = med_embedding_seg_precision

    medians = [md_encoding_true_positives, md_embedding_true_positives,
               md_encoding_precisions, md_embedding_precisions]
    return medians

def compute_average_segment_statistic(args, query_IDs, database_IDs):
    # average
    av_encoding_true_positives = dict()
    av_embedding_true_positives = dict()
    av_encoding_precisions = dict()
    av_embedding_precisions = dict()

    for s in range(int(args.num_seg)):
        print('computing metrics for segment ...' + str(s))

        # get distances for top k
        for segment_file in os.listdir(args.data_dir):
            f = os.path.join(args.data_dir, segment_file)
            base = 'seg.' + str(s)
            # segment's plink file
            if (base in f and args.plink_ext in f):
                s_plink_file = f
                # plink dictionary
                plink_dict = basic_datastructures.get_plink_distances(
                    database_IDs, query_IDs, s_plink_file)
                plink_top = compute_conditions.compute_top_k(plink_dict, int(args.k))

            # segment's faiss file for embeddings
            elif (base in f and args.faiss_enc_ext in f):
                s_faiss_enc_file = f
                encoding_faiss_dict = basic_datastructures.get_faiss_distances(
                    database_IDs, query_IDs, s_faiss_enc_file)
                encoding_faiss_top = compute_conditions.compute_top_k(encoding_faiss_dict, int(args.k))

            # segment's faiss file for encodings
            elif (base in f and args.faiss_emb_ext in f):
                s_faiss_emb_file = f
                embedding_faiss_dict = basic_datastructures.get_faiss_distances(
                    database_IDs, query_IDs, s_faiss_emb_file)
                embedding_faiss_top = compute_conditions.compute_top_k(embedding_faiss_dict, int(args.k))

        # compute metrics
        # get all true positives
        encoding_true_positives = compute_conditions.compute_true_positives(plink_top, encoding_faiss_top)
        embedding_true_positives = compute_conditions.compute_true_positives(plink_top, embedding_faiss_top)
        encoding_seg_TP_values = [encoding_true_positives[p] for p in query_IDs]
        embedding_seg_TP_values = [embedding_true_positives[p] for p in query_IDs]
        # get avg of true positives
        avg_encoding_seg_TP = (sum(encoding_seg_TP_values) / len(encoding_seg_TP_values))
        avg_embedding_seg_TP = (sum(embedding_seg_TP_values) / len(embedding_seg_TP_values))

        # get all precision
        encoding_precisions = compute_conditions.confusion_matrix(encoding_true_positives, int(args.k))
        embedding_precisions = compute_conditions.confusion_matrix(embedding_true_positives, int(args.k))
        encoding_seg_prec_values = [encoding_precisions[p] for p in query_IDs]
        embedding_seg_prec_values = [embedding_precisions[p] for p in query_IDs]
        # get avg of precision
        avg_encoding_seg_precision = (sum(encoding_seg_prec_values) / len(encoding_seg_prec_values))
        avg_embedding_seg_precision = (sum(embedding_seg_prec_values) / len(embedding_seg_prec_values))

        # add averages to dictionary
        av_encoding_true_positives[s] = avg_encoding_seg_TP
        av_embedding_true_positives[s] = avg_embedding_seg_TP
        av_encoding_precisions[s] = avg_encoding_seg_precision
        av_embedding_precisions[s] = avg_embedding_seg_precision

    averages = [av_encoding_true_positives, av_embedding_true_positives,
                av_encoding_precisions, av_embedding_precisions]
    return averages


if __name__ == '__main__':
    main()
