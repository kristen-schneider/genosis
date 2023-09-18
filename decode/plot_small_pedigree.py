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
    ss_sample_results_dir = args.ss_sample_results_dir + '/'
    out_dir = args.out_dir + '/'

    sample_IDs = sp.read_sample_IDs(sample_IDs_txt)
    for sample_ID in sample_IDs:
        print(sample_ID)
        r_file = ss_sample_results_dir + sample_ID + '.csv'
        results = read_results_file(out_dir, r_file)

def read_results_file(out_dir, r_file):
    results = dict()
    r_file = open(r_file, 'r')
    sample_ID = r_file.readline().strip().split('_')[0]
    x_labels = []
    y_scores = []

    for line in r_file:
        if '_' in line:
            pass
        else:
            match_ID = line.strip().split(':')[0]
            try:
                results[match_ID] += float(line.strip().split(':')[1])
            except KeyError:
                results[match_ID] = float(line.strip().split(':')[1])

    for match_ID in results:
        x_labels.append(match_ID)
        y_scores.append(results[match_ID])
        
    plt.figure(figsize=(20, 10))
    plt.bar(x_labels, y_scores, color='olivedrab')
    plt.title('Query ID: ' + sample_ID)
    plt.xticks(rotation=90)
    plt.xlabel('Match IDs')
    plt.ylabel('Scores')
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.savefig(out_dir + sample_ID + '.png')
    plt.close()

if __name__ == '__main__':
    main()
