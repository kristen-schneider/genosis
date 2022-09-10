import argparse
import matplotlib.pyplot as plt
import sys, os

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--faiss_encoding')
    parser.add_argument('--faiss_embedding')
    parser.add_argument('--plink')
    parser.add_argument('--png')

    return parser.parse_args()

def main():
    args = get_args()

    # plot averages
    faiss_enc = read_data_file(args.faiss_encoding)
    faiss_emb = read_data_file(args.faiss_embedding)
    plink = read_data_file(args.plink)

    plot_segments(faiss_enc, faiss_emb, plink, args.png)

def plot_segments(faiss_enc, faiss_emb, plink, png_file):
    num_segments = len(faiss_enc)
    x = range(num_segments)

    faiss_enc_y = [faiss_enc[v] for v in faiss_enc.keys()]
    faiss_emb_y = [faiss_emb[v] for v in faiss_emb.keys()]
    plink_y = [plink[v] for v in plink.keys()]

    plt.figure(figsize=(20, 10))
    ax1 = plt.subplot(111)
    # True Positive
    tp_title_str = 'Query Times for ' + str(num_segments) + ' Segments'
    ax1.set_title(tp_title_str, fontsize=30)
    faiss_enc_plot = ax1.plot(x, faiss_enc_y, color='olivedrab', linestyle='-', marker='o')
    faiss_emb_plot = ax1.plot(x, faiss_emb_y, color='firebrick', linestyle='-', marker='x')
    plink_plot = ax1.plot(x, plink_y, color='steelblue', linestyle='-', marker='x')

    # ax1.set_xticks(range(num_segments))
    ax1.set_xlabel('Segment', fontsize=20)
    ax1.set_ylabel('Time (seconds)', fontsize=20)
    ax1.spines.right.set_visible(False)
    ax1.spines.top.set_visible(False)

    plt.legend(['FAISS encodings', 'FIASS embeddings', 'PLINK'], prop={'size': 30})
    plt.savefig(png_file)


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

