import argparse
import plot_decode as pd

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
        super_pops, sub_pops = 
        super_labels, sub_labels = 
        
        super_ancestry, sub_ancestry = label_acestry(ancestry_file)
        


def get_relation_scores(samples, pair_relations, data_dir, ext):
   
    all_relations_scores = {}

    for sample in samples:
        #print(sample)
        # open sample file
        sample_file = data_dir + '/' + sample + ext + '.csv'
        query = sample
        with open(sample_file, 'r') as file:
            # combine any sample amtches across haps
            sample_haps = {}
            for line in file:
                if ":" not in line:
                    hap = line
                else:
                    values = line.strip().split()
                    match_ID = values[0].split(":")[0][:-2]
                    score = float(values[1])
                    try:
                        sample_haps[match_ID].append(score)
                    except KeyError:
                        sample_haps[match_ID] = [score]
        
        # for all matches with this sample append and get score
        for m in sample_haps:
            curr_score = sum(sample_haps[m])
            relation = pair_relations[query][m]
            if relation == "sibling" and ext == '.pop':
                if curr_score > 100: print(query, m, sample_haps[m], curr_score)
            try:
                all_relations_scores[relation].append(curr_score)
            except KeyError:
                all_relations_scores[relation] = [curr_score]
    return all_relations_scores
            
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
    
def label_ancestry(ancestry_file):
    super_ancestry = {}
    sub_ancestry = {}

    with open(ancestry_file, 'r') as file:
        for line in file:
            values = line.strip().split(',')
            sample_ID = values[0]
            sub_pop = values[3]
            super_pop = values[5]
                   
            super_ancestry[sample_ID] = super_pop
            sub_ancestry[sample_ID] = sub_pop

    return super_ancestry, sub_ancestry
    
def write_relations_scores(all_relations_scores, out_f):
    o = open(out_f, 'w')
    for relation in all_relations_scores:
        o.write(relation + ',' + str(all_relations_scores[relation]))
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
