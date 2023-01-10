import argparse

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--sample_info')
    parser.add_argument('--pop')
    parser.add_argument('--out_dir')
    return parser.parse_args()

def main():
    args = get_args()
    sample_info = args.sample_info
    pop = args.pop
    out_dir = args.out_dir

    get_subpopulation(sample_info, pop, out_dir)

def get_subpopulation(sample_info, pop, out_dir):

    in_stream = open(sample_info, 'r')
    out_stream = open(out_dir + pop + '_samples.txt', 'w')
    
    for line in in_stream:
        L = line.strip().split('\t')
        population_code = L[3]
        superpopulation_code = L[5]
        if superpopulation_code == pop:
            out_stream.write(line)
            #print(population_code, superpopulation_code)

    in_stream.close()
    out_stream.close()




if __name__ == '__main__':
    main()
