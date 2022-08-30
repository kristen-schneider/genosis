import argparse
import matplotlib.pyplot as plt
import sys, os

import os, sys, inspect

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir)

import compute_conditions
import basic_datastructures


'''
yeah, I think it'd just be easier to
count the TP/FP/FN and do average/stddev
of prec/recall across all queries.
That would make it easier to quatifiably
compare with the encoding based performance.
'''

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--data_dir')
    parser.add_argument('--plink')
    parser.add_argument('--faiss_enc')
    parser.add_argument('--faiss_emb')
    parser.add_argument('--k')
    parser.add_argument('--train')
    parser.add_argument('--test')
    parser.add_argument('--outdir')
    parser.add_argument('--num_seg')
    return parser.parse_args()

def main():
    args = get_args()
    query_IDs = basic_datastructures.get_db_q_IDs(args.test)
    database_IDs = basic_datastructures.get_db_q_IDs(args.train)

    [ENC_TP, EMB_TP, ENC_P, EMB_P] = many_segments(args, query_IDs, database_IDs, args.num_seg)
    seg_precision_png = args.outdir+'many_seg_precision.png'
    plot_many_segments(ENC_P, EMB_P, seg_precision_png, 'precision')
    return 0

def many_segments(args, query_IDs, database_IDs, number_segments):
    av_encoding_true_positives = dict()
    av_embedding_true_positives = dict()
    av_encoding_precisions = dict()
    av_embedding_precisions = dict()

    for s in range(int(number_segments)):
        encoding_segment_TP = -1
        embedding_segment_TP = -1
        encoding_segment_precision = -1
        embedding_segment_precision = -1

        for segment_file in os.listdir(args.data_dir):
            f = os.path.join(args.data_dir, segment_file)
            base = 'seg.'+ str(s)
            # segment's plink file
            if ( base in f and args.plink in f):
                s_plink_file = f
                # plink dictionary
                plink_dict = basic_datastructures.get_plink_distances(
                    database_IDs, query_IDs, s_plink_file)
                plink_top = compute_conditions.compute_top_k(plink_dict, int(args.k))

            # segment's faiss file for embeddings
            elif( base in f and args.faiss_enc in f):
                s_faiss_enc_file = f
                encoding_faiss_dict = basic_datastructures.get_faiss_distances(
                    database_IDs, query_IDs, s_faiss_enc_file)
                encoding_faiss_top = compute_conditions.compute_top_k(encoding_faiss_dict, int(args.k))

            # segment's faiss file for encodings
            elif( base in f and args.faiss_emb in f):
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
        # get sum of true positives
        encoding_segment_TP = (sum(encoding_seg_TP_values)/len(encoding_seg_TP_values))
        embedding_segment_TP = (sum(embedding_seg_TP_values)/len(embedding_seg_TP_values))

        encoding_precisions = confusion_matrix(encoding_true_positives, int(args.k))
        embedding_precisions = confusion_matrix(embedding_true_positives, int(args.k))
        encoding_seg_prec_values = [encoding_precisions[p] for p in query_IDs]
        embedding_seg_prec_values = [embedding_precisions[p] for p in query_IDs]
        # get sum of true positives
        encoding_segment_precision = (sum(encoding_seg_prec_values)/len(encoding_seg_prec_values))
        embedding_segment_precision = (sum(embedding_seg_prec_values)/len(embedding_seg_prec_values))

        av_encoding_true_positives[s] = encoding_segment_TP
        # except KeyError: av_encoding_true_positives[s] = [encoding_segment_TP]

        av_embedding_true_positives[s] = embedding_segment_TP
        # except KeyError: av_embedding_true_positives[s] = [embedding_segment_TP]

        av_encoding_precisions[s] = encoding_segment_precision
        # except KeyError: av_encoding_precisions[s] = [encoding_segment_precision]

        av_embedding_precisions[s] = embedding_segment_precision
        # except KeyError: av_embedding_precisions[s] = embedding_segment_precision]

    return [av_encoding_true_positives, av_embedding_true_positives,
            av_encoding_precisions, av_embedding_precisions]



def single_segment(query_IDs, database_IDs, args):
    # plink dictionary
    plink_dict = basic_datastructures.get_plink_distances(database_IDs, query_IDs, args.plink)
    encoding_faiss_dict = basic_datastructures.get_faiss_distances(
        database_IDs, query_IDs, args.faiss_enc)
    embedding_faiss_dict = basic_datastructures.get_faiss_distances(
        database_IDs, query_IDs, args.faiss_emb)

    plink_top = compute_conditions.compute_top_k(plink_dict, int(args.k))
    encoding_faiss_top = compute_conditions.compute_top_k(encoding_faiss_dict, int(args.k))
    embedding_faiss_top = compute_conditions.compute_top_k(embedding_faiss_dict, int(args.k))


    # positives = compute_conditions.compute_positives(plink_top, faiss_top)
    encoding_true_positives = compute_conditions.compute_true_positives(plink_top, encoding_faiss_top)
    embedding_true_positives = compute_conditions.compute_true_positives(plink_top, embedding_faiss_top)

    encoding_query_precisions = confusion_matrix(encoding_true_positives, int(args.k))
    embedding_query_precisions = confusion_matrix(embedding_true_positives, int(args.k))

    plot_hist(encoding_true_positives, embedding_true_positives,
              args.outdir+'TP-hist.png', 'True Positive')
    plot_hist(encoding_query_precisions, embedding_query_precisions,
              args.outdir+'precision-hist.png', 'Precision')

def confusion_matrix(TP, k):
    #    ___P________N___
    # P |___TP___|___FN___|
    # N |___FP___|___TN___|
    query_precisions = dict()

    for query_ID in TP.keys():
        q_TP = TP[query_ID]
        q_FP = k - q_TP
        q_FN = q_FP
        # q_TN = NA

        # precision = TP/(TP+FP) = recall
        q_precision = q_TP/(q_TP+q_FP)

        query_precisions[query_ID] = q_precision

    return query_precisions

def plot_many_segments(seg_av_encoding, seg_av_embedding, out_png, metric):

    x = range(len(seg_av_encoding))
    enc_y = [seg_av_encoding[v] for v in seg_av_encoding.keys()]
    emb_y = [seg_av_embedding[v] for v in seg_av_embedding.keys()]

    plt.figure(figsize=(30, 20))
    ax1 = plt.subplot(111)
    # title_str = metric + ' for ' + str(len(query_IDs)) + ' Queries'
    # ax1.set_title(title_str, fontsize=30)
    ax1.plot(x, enc_y, color='olivedrab', linestyle='-', marker='o')
    ax1.plot(x, emb_y, color='salmon', linestyle='-', marker='x')
    ax1.set_xlabel('Segment', fontsize=20)
    ax1.set_ylabel(metric, fontsize=20)
    ax1.spines.right.set_visible(False)
    ax1.spines.top.set_visible(False)
    plt.savefig(out_png)


def plot_hist(enc_query_precision, emb_query_precision, out_png, metric):
    # Positive Predictive Value (PPV)
    # Precision
    # (TP / (TP + FP))
    query_IDs = enc_query_precision.keys()
    x = range(len(query_IDs))
    encoding_precision = [enc_query_precision[p] for p in query_IDs]
    embedding_precision = [emb_query_precision[p] for p in query_IDs]

    plt.figure(figsize=(30, 20))
    ax1 = plt.subplot(211)
    title_str = metric+' for '+str(len(query_IDs))+' Queries'
    ax1.set_title(title_str, fontsize=30)
    ax1.hist(encoding_precision, color='olivedrab', alpha=0.8)
    ax1.set_xlabel('Encoding '+metric, fontsize=20)
    ax1.set_ylabel('Frequency', fontsize=20)
    ax1.spines.right.set_visible(False)
    ax1.spines.top.set_visible(False)

    ax2 = plt.subplot(212)
    ax2.hist(embedding_precision, color='olivedrab', alpha=0.8)
    ax2.set_xlabel('Embedding '+metric, fontsize=20)
    ax2.set_ylabel('Frequency', fontsize=20)
    ax2.spines.right.set_visible(False)
    ax2.spines.top.set_visible(False)

    plt.savefig(out_png)


def plot_precision_bar(enc_query_precision, emb_query_precision, out_png):
    # Positive Predictive Value (PPV)
    # Precision
    # (TP / (TP + FP))
    query_IDs = enc_query_precision.keys()
    encoding_precision = [enc_query_precision[p] for p in query_IDs]
    embedding_precision = [emb_query_precision[p] for p in query_IDs]

    plt.figure(figsize=(15, 15))
    ax1 = plt.subplot(111)
    y = [(sum(encoding_precision)/len(encoding_precision)),
         (sum(embedding_precision)/len(embedding_precision))]
    x = range(len(y))

    ax1.bar(x, y, color=['yellowgreen', 'salmon'])

    ax1.set_xticks(x, ['encoding', 'embedding'], rotation=0, fontsize=20)
    # ax1.set_xlabel('Query ID', fontsize=20)
    ax1.set_ylabel('Precision', fontsize=20)
    ax1.spines.right.set_visible(False)
    ax1.spines.top.set_visible(False)

    plt.title('Average Precision for all queries '
              '\n one seg = 5000snps', fontsize=30)
    plt.savefig(out_png)

def plot_precision(enc_query_precision, emb_query_precision, out_png):
    # Positive Predictive Value (PPV)
    # Precision
    # (TP / (TP + FP))
    query_IDs = enc_query_precision.keys()
    x = range(len(query_IDs))
    encoding_precision = [enc_query_precision[p] for p in query_IDs]
    embedding_precision = [emb_query_precision[p] for p in query_IDs]

    plt.figure(figsize=(40, 15))
    ax1 = plt.subplot(111)
    #
    ax1.scatter(x, encoding_precision, color='yellowgreen', marker='o')
    ax1.scatter(x, embedding_precision, color='salmon', marker='x')

    ax1.set_xticks(x, query_IDs, rotation=90)
    ax1.set_xlabel('Query ID', fontsize=20)
    ax1.set_ylabel('Precision', fontsize=20)
    ax1.spines.left.set_linewidth(4)
    ax1.spines.right.set_visible(False)
    ax1.spines.top.set_visible(False)

    plt.savefig(out_png)



# def plot_ROC():
#     # x = False Positive Rate (FP / N) = 1 - TNR
#     # y = True Positive Rate (TPR) = (TP / P) = 1 - FNR
#
# def plot_accuracy():
#     # Accuracy (ACC).
#     # (TP + TN) / (P + N)
#
# def plot_TPR():
#     # True Positive Rate (TPR).
#     # Recall. Sensitivity.
#     # (TP / P) = 1 - FNR
# def plot


if __name__ == '__main__':
    main()
