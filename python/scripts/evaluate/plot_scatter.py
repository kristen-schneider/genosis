import argparse
from collections import defaultdict
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch


def parse_args():
    parser = argparse.ArgumentParser(description='scatter plot ibd vs our score')
    parser.add_argument('--ibd', type=str, help='ground truth IBD file')
    parser.add_argument('--knn', type=str, help='directory with knn scores')
    parser.add_argument('--samples', type=str, help='file with sample names')
    parser.add_argument('--relations', type=str, help='file with named relations')
    parser.add_argument('--png', type=str, help='out png file')
    return parser.parse_args()

def main():
    args = parse_args()

    ibd_file = args.ibd
    knn_dir = args.knn
    samples_file = args.samples
    relations_file = args.relations
    png = args.png

    samples = read_samples(samples_file)
    relations = read_relations(relations_file)
    svs_scores = knn_pairwise_scores(samples, knn_dir, '.pop')
    ibd_scores = ibd_pairwise_scores(ibd_file)
    print('plotting')
    scatter_plot(samples, svs_scores, ibd_scores, relations, png)

def read_relations(relations_file):
    relations = {}
    with open(relations_file, 'r') as file:
        for line in file:
            values = line.strip().split(',')
            familyID = values[0]
            i1 = values[1]
            i2 = values[2]
            relation = values[3]
            try:
                relations[i1].update({i2: relation})
            except KeyError:
                relations[i1] = {i2: relation}
    return relations    


def scatter_plot(samples, svs_scores, ibd_scores, relations, png):
    svs_plot_scores = []
    ibd_plot_scores = []
    colors = []

    colors_key = {'self':'black',
                    'child': 'olivedrab', 'parent': 'olivedrab',
                    'sibling': 'red',
                    'grandchild': 'steelblue', 'grandparent': 'steelblue',
                    'niece-nephew': 'gold', 'aunt-uncle': 'gold',
                    'great-grandchild': 'lightblue', 'great-grandparent': 'lightblue',
                    '1-cousin': 'hotpink',
                    'great-niece-nephew': 'khaki', 'great-aunt-uncle': 'khaki',
                    '1-cousin-1-removed': 'pink',
                    '2-cousin': 'lightpink',
                    'unrelated': 'gray'}
    legend_elements = []
    for c in colors_key:
        legend_elements.append(Line2D([0], [0], marker='o', color=colors_key[c], label=c,
                          markerfacecolor=colors_key[c], markersize=15))
    
    for q in samples:
        for m in samples:
            try:
                ibd_score = ibd_scores[q][m]
                try:
                    svs_score = svs_scores[q][m]
                    try:
                        ibd_plot_scores.append(ibd_scores[q][m])
                        svs_plot_scores.append(svs_scores[q][m])
                        relation = relations[q][m]
                        try:
                            colors.append(colors_key[relation])
                        except KeyError:
                            colors.append('black')
                    except KeyError:
                        pass
                except KeyError:
                    pass
            except KeyError:
                pass
                
    print(len(svs_plot_scores), len(ibd_plot_scores))

    fig, ax = plt.subplots(figsize=(20,15))

    ax.scatter(svs_plot_scores, ibd_plot_scores, c=colors)
    ax.set_xlabel("svs scores")
    ax.set_ylabel("ibd_scores")
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    

    ax.legend(handles=legend_elements)
    plt.savefig(png)

    
def ibd_pairwise_scores(ibd_file):
   
    all_pairwise_scores = defaultdict(dict)

    with open(ibd_file, 'r') as file:
        # combine any sample matches across haps
        for line in file:
            values = line.strip().split()
            query_ID = values[0]
            match_ID = values[1]
            score = float(values[8])
            
            try:
                all_pairwise_scores[query_ID][match_ID] += score
            except KeyError:
                all_pairwise_scores[query_ID][match_ID] = score
        
    return all_pairwise_scores
 

def knn_pairwise_scores(samples, data_dir, ext):
   
    all_pairwise_scores = defaultdict(dict)

    for sample in samples:
        # open sample file
        sample_file = data_dir + sample + ext + '.csv'
        query_ID = sample
        with open(sample_file, 'r') as file:
            # combine any sample matches across haps
            for line in file:
                if ":" not in line:
                    hap = line
                else:
                    values = line.strip().split()
                    match_ID = values[0].split(":")[0][:-2]
                    score = float(values[1])
                    try:
                        all_pairwise_scores[query_ID][match_ID] += score
                    except KeyError:
                        all_pairwise_scores[query_ID][match_ID] = score
        
    return all_pairwise_scores

def read_samples(samples_file):
    samples = []
    with open(samples_file, 'r') as file:
        for line in file:
            samples.append(line.strip())
    return samples


if __name__ == '__main__':
    main()
