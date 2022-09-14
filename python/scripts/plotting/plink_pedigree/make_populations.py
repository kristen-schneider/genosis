import argparse

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--sample_data')
    parser.add_argument('--samples_dir')
    return parser.parse_args()

def main():
    args = get_args()
    pop_samples = make_pop_samples_dict(args.sample_data)
    pop_codes = pop_samples.keys()

    for pc in pop_codes:
        o = open(args.samples_dir + pc + '.samples', 'w')
        o.write(pc + '\n')
        for s in pop_samples[pc]:
            o.write(s + '\n')
        o.close()



def make_pop_samples_dict(sample_metadata): 
    pop_samples = dict()
    
    f = open(sample_metadata, 'r')
    header = f.readline()
    for line in f:
        L = line.strip().split()
        sample_ID = L[0]
        pop_code = L[3]
        try:
            pop_samples[pop_code].append(sample_ID)
        except KeyError:
            pop_samples[pop_code] = [sample_ID]

    f.close()
    return pop_samples

if __name__ == '__main__':
    main()
