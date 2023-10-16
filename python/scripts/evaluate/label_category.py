import argparse
import plot_violin as vp
import plot_scatter as sp
import pedigree
import ancestry

def parse_args():
    parser = argparse.ArgumentParser(description='Label pedigree relations between individuals')
    parser.add_argument('--relations', type=str, help='relations file')
    parser.add_argument('--ancestry', type=str, help='ancestry file')
    parser.add_argument('--samples', type=str, help='samples_file')
    parser.add_argument('--data_dir', type=str, help='data directory for out csv files')
    return parser.parse_args()

def main():
    args = parse_args()
    relations_file = args.relations
    ancestry_file = args.ancestry
    samples_file = args.samples
    data_dir = args.data_dir
   
     
    samples = read_samples(samples_file)
    
    # if peidgree --> relations file
    if ancestry_file == 'None':
        print("evaluating pedigree data...")

        pair_relations = read_relations(relations_file)

        all_relations_scores_svs = get_relation_scores(samples, pair_relations, data_dir, ext='.svs')
        all_relations_scores_pop = get_relation_scores(samples, pair_relations, data_dir, ext='.pop')
        all_relations_scores_ibd = get_relation_scores(samples, pair_relations, data_dir, ext='.ibd')

        write_relations_scores(all_relations_scores_svs, data_dir + '/SVS.results')    
        write_relations_scores(all_relations_scores_pop, data_dir + '/POP.results')    
        write_relations_scores(all_relations_scores_ibd, data_dir + '/IBD.results')    

        pd.violin_plot(all_relations_scores_svs, "SVS scores", data_dir + '/SVS.png')
        pd.violin_plot(all_relations_scores_pop, "popcount", data_dir + '/POP.png')
        pd.violin_plot(all_relations_scores_ibd, "merge-gap", data_dir + '/IBD.png')
	
    # else --> ancestry file    
    else:
        print("evaluating ancestry data...")
        super_pops_sub_pops = get_super_sub(ancestry_file)
        super_labels, sub_labels = get_ancestry_names(ancestry_file)
        
        super_ancestry, sub_ancestry = label_acestry(ancestry_file)
        
            
def write_category_scores(all_category_scores, out_f):
    o = open(out_f, 'w')
    for category in all_category_scores:
        o.write(category + ',' + str(all_category_scores[category]))
        o.write('\n')
    o.close()

def read_samples(samples_file):
    samples = []
    with open(samples_file, 'r') as file:
        for line in file:
            samples.append(line.strip())
    return samples

if __name__ == '__main__':
    main()
