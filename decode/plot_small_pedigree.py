import argparse
from collections import defaultdict
import matplotlib.pyplot as plt
import small_pedigree as sp

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--sample_IDs')
    parser.add_argument('--ss_sample_results_dir')
    parser.add_argument('--out_dir')

    return parser.parse_args()

def main():
    args = get_args()
    sample_IDs_txt = args.sample_IDs
    ss_sample_results_dir = args.ss_sample_results_dir
    out_dir = args.out_dir

    sample_IDs = sp.read_sample_IDs(sample_IDs_txt)
    for sample_ID in sample_IDs:
        print(sample_ID)
        r_file = ss_sample_results_dir + sample_ID + '.csv'
        results = read_results_file(out_dir, r_file)

def read_results_file(out_dir, r_file):
    results = defaultdict(dict)
    r_file = open(r_file, 'r')
    hap0 = ''
    hap1 = ''
    x_labels = []
    y_scores = []

    for line in r_file:
        if '_' in line:
            if hap0 == '':
                hap0 = line.strip()
            else:
                # plot a bar graph of the scores for each match and the x labels are the match IDs
                plt.figure(figsize=(20, 10))
                plt.bar(x_labels, y_scores, color='olivedrab')
                plt.title('Query ID: ' + hap0)
                plt.xticks(rotation=90)
                plt.xlabel('Match IDs')
                plt.ylabel('Scores')
                plt.gca().spines['top'].set_visible(False)
                plt.gca().spines['right'].set_visible(False)
                plt.savefig(out_dir + hap0 + '.png')
                plt.close()
                x_labels = []
                y_scores = []
                hap1 = line.strip()

        else:
            x_labels.append(line.strip().split(':')[0])
            y_scores.append(line.strip().split(':')[1])

    plt.figure(figsize=(20, 10))
    plt.bar(x_labels, y_scores, color='olivedrab')
    plt.title('Query ID: ' + hap0)
    plt.xticks(rotation=90)
    plt.xlabel('Match IDs')
    plt.ylabel('Scores')
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.savefig(out_dir + hap1 + '.png')
    plt.close()
    x_labels = []
    y_scores = []
    hap1 = line.strip()

if __name__ == '__main__':
    main()
