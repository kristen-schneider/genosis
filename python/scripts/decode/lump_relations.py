import argparse
import plot_decode as pd

def parse_args():
    parser = argparse.ArgumentParser(description='Label pedigree relations between individuals')
    parser.add_argument('--relations', type=str, help='relations file')
    parser.add_argument('--samples', type=str, help='samples_file')
    parser.add_argument('--data_dir', type=str, help='data directory for out csv files')
    return parser.parse_args()

def main():
    args = parse_args()
    relations_file = args.relations
    samples_file = args.samples
    data_dir = args.data_dir
    
    samples = read_samples(samples_file)
    pair_relations = read_relations(relations_file)

    all_relations_scores_pop = get_relation_scores(samples, pair_relations, data_dir, ext='.pop')
    all_relations_scores_ibd = get_relation_scores(samples, pair_relations, data_dir, ext='.ibd')

    write_relations_scores(all_relations_scores_pop, data_dir + '/POP.results')    
    write_relations_scores(all_relations_scores_pop, data_dir + '/IBD.results')    

    pd.violin_plot(all_relations_scores_pop, data_dir + '/POP.png')
    pd.violin_plot(all_relations_scores_ibd, data_dir + '/IBD.png')
	
def get_relation_scores(samples, pair_relations, data_dir, ext):
   
    all_relations_scores = {}

    for sample in samples:
        print(sample)
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
                    match_ID = values[0].split(":")[0]
                    score = float(values[1])
                    try:
                        sample_haps[match_ID] += score
                    except KeyError:
                        sample_haps[match_ID] = score
        # for all matches with this sample append and get score
        for m in sample_haps:
            curr_score = sample_haps[m]
            relation = pair_relations[query][m]
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
            i1 = values[0]
            i2 = values[1]
            relation = values[2]
            try:
		relations[i1].update({i2: relation})
            except KeyError:
                relations[i1] = {i2: relation}
    return relations
    
def read_samples(samples_file):
    samples = []
    with open(samples_file, 'r') as file:
        for line in file:
            samples.append(line.strip())
    return samples
    
def write_relations_scores(all_relations_scores, out_f):
    o = open(out_f, 'w')
    for relation in all_relations_scores:
        o.write(relation + ',' + str(all_relations_scores[relation]))
	o.write('\n')
   o.close()


if __name__ == '__main__':
	    main()
