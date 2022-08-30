import argparse
import matplotlib.pyplot as plt
import sys, os

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--avg_tp_encoding')
    parser.add_argument('--avg_tp_embedding')
    parser.add_argument('--avg_p_encoding')
    parser.add_argument('--avg_p_embedding')
    parser.add_argument('--med_tp_encoding')
    parser.add_argument('--med_tp_embedding')
    parser.add_argument('--med_p_encoding')
    parser.add_argument('--med_p_embedding')
    parser.add_argument('--avg_png')
    parser.add_argument('--med_png')

    return parser.parse_args()

def main():
    args = get_args()

    # plot averages
    avg_tp_enc = read_data_file(args.avg_tp_encoding)
    avg_tp_emb = read_data_file(args.avg_tp_embedding)
    avg_p_enc = read_data_file(args.avg_p_encoding)
    avg_p_emb = read_data_file(args.avg_p_embedding)
    plot_segments(avg_tp_enc, avg_tp_emb, avg_p_enc, avg_p_emb, 'Average', args.avg_png)

    # plot medians
    med_tp_enc = read_data_file(args.med_tp_encoding)
    med_tp_emb = read_data_file(args.med_tp_embedding)
    med_p_enc = read_data_file(args.med_p_encoding)
    med_p_emb = read_data_file(args.med_p_embedding)
    plot_segments(med_tp_enc, med_tp_emb, med_p_enc, med_p_emb, 'Median', args.med_png)


def plot_segments(tp_encoding, tp_embedding, p_encoding, p_embedding, stat, outpng):

    num_segments = len(tp_encoding)
    x = range(num_segments)

    tp_enc_y = [tp_encoding[v] for v in tp_encoding.keys()]
    tp_emb_y = [tp_embedding[v] for v in tp_embedding.keys()]
    p_enc_y = [p_encoding[v] for v in p_encoding.keys()]
    p_emb_y = [p_embedding[v] for v in p_embedding.keys()]

    plt.figure(figsize=(30, 20))
    ax1 = plt.subplot(211)
    # True Positive
    tp_title_str = stat + ' True Positive and Precision for ' + str(num_segments) + ' Segments'
    ax1.set_title(tp_title_str, fontsize=30)
    tp_enc_plot = ax1.plot(x, tp_enc_y, color='olivedrab', linestyle='-', marker='o')
    tp_emb_plot = ax1.plot(x, tp_emb_y, color='salmon', linestyle='-', marker='x')
    ax1.set_xticks([])
    # ax1.set_xlabel('Segment', fontsize=20)
    ax1.set_ylabel('True Positive', fontsize=20)
    ax1.spines.right.set_visible(False)
    ax1.spines.top.set_visible(False)
    ax1.spines.bottom.set_visible(False)


    # Precision
    ax2 = plt.subplot(212)
    p_enc_plot = ax2.plot(x, p_enc_y, color='olivedrab', linestyle='-', marker='o')
    p_emb_plot = ax2.plot(x, p_emb_y, color='salmon', linestyle='-', marker='x')
    ax2.set_xticks(range(num_segments), fontsize=15)
    ax2.set_xlabel('Segment', fontsize=20)
    ax2.set_ylabel('Precision', fontsize=20)
    ax2.spines.right.set_visible(False)
    ax2.spines.top.set_visible(False)

    plt.legend(['ENCODINGS', 'EMBEDDINGS'], prop={'size': 30})
    plt.savefig(outpng)

def read_data_file(in_file):
    data_dict = dict()
    f = open(in_file, 'r')
    for line in f:
        l = line.strip().split()
        segment_idx = int(l[0])
        metric = float(l[1])
        data_dict[segment_idx] = metric
    f.close()
    return data_dict

if __name__ == '__main__':
    main()
