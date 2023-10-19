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

        # create a dictionary of samples pairs and their pedigree relationship label
        pair_relations = read_relations(relations_file)
        
        # create dictionaries of scores for different peigree relationship labels
        all_pedigree_scores_svs = pedigree.get_pedigree_scores(samples, pair_relations, data_dir, ext='.svs')
        all_pedigree_scores_pop = pedigree.get_pedigree_scores(samples, pair_relations, data_dir, ext='.pop')
        all_pedigree_scores_ibd = pedigree.get_pedigree_scores(samples, pair_relations, data_dir, ext='.ibd')

        # write these scores to output files for saving/transfering
        write_category_scores(all_pedigree_scores_svs, data_dir + '/SVS.results')    
        write_category_scores(all_pedigree_scores_pop, data_dir + '/POP.results')    
        write_category_scores(all_relations_scores_ibd, data_dir + '/IBD.results')    
    
        # violin plot of scores
        pv.violin_plot(all_pedigree_scores_svs, "SVS scores", data_dir + '/SVS.png')
        pv.violin_plot(all_pedigree_scores_pop, "popcount", data_dir + '/POP.png')
        pv.violin_plot(all_pedigree_scores_ibd, "merge-gap", data_dir + '/IBD.png')
	
    # else --> ancestry file    
    else:
        print("evaluating ancestry data...")
        
        # dictionoary of superpop1: [subpop1, subpop2, ...]
        super_pops_sub_pops = ancestry.get_super_sub(ancestry_file)
        # dictionaries of ancestryID: ancestry name
        super_labels, sub_labels = ancestry.get_ancestry_names(ancestry_file)
        # dictionaries of samppleID: ancestryID
        super_ancestry, sub_ancestry = ancestry.label_acestry(ancestry_file)
            
        # get ancestry scores
            
            

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
