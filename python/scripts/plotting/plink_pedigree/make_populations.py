import argparse

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--sample_data')
    parser.add_argument('--samples_dir')
    return parser.parse_args()

def main():
    args = get_args()
    make_pop_summary_file(args.sample_data, args.samples_dir)
    #pop_samples = make_pop_samples_dict(args.sample_data)

def make_pop_summary_file(sample_metadata, out_dir):
    o = open(out_dir + 'population.summary', 'w')

    f = open(sample_metadata, 'r')
    header = f.readline()
    for line in f:
        L = line.strip().split('\t')
        sample_ID = L[0]
        pop_code = L[3]
        pop_name = L[4]
        super_pop_code = L[5]
        super_pop_name = L[6]
        
        out_line = sample_ID + ',' + pop_code + ',' + pop_name + ',' + super_pop_code + ',' + super_pop_name + '\n'
        o.write(out_line)
    f.close()
    o.close()

def write_pop_samples(pop_samples):
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
