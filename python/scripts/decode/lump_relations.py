import argparse
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
    
    all_relations_scores = {}
    
    samples = read_samples(samples_file)
    pair_relations = read_relations(relations_file)

def get_relation_scores(samples, pair_relations, data_dir, all_relations_scores):
    for sample in samples:
        print(sample)
        # open sample file
        sample_file = data_dir + '/' + sample + '.csv'
        query = sample
        with open(sample_file, 'r') as file:
            for line in file:
                if ":" not in line:
                    break
                values = line.strip().split()
                match = values[0].split(":")[0]
                score = float(values[1])
                relation = pair_relations[query][match]
                try:
                    all_relations_scores[relation].append(score)
                except KeyError:
                    all_relations_scores[relation] = [score]
            
def read_relations(relations_file):
    relations = {}
    with open(relations_file, 'r') as file:
        for line in file:
            values = line.strip().split()
            i1 = values[0]
            i2 = values[1]
            relation = values[2]
            relations[i1] = {i2: relation}
    return relations
    
def read_samples(samples_file):
    samples = []
    with open(samples_file, 'r') as file:
        for line in file:
            samples.append(line.strip())
    return samples
    
if __name__ == '__main__':
    main()
