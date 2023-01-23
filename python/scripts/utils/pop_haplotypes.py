import argparse

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--sample_IDs')
    parser.add_argument('--all_haplotypes')
    return parser.parse_args()

def main():
    args = get_args()
    sample_IDs = args.sample_IDs
    all_haplotypes = args.all_haplotypes

    get_subpopulation_haps(sample_IDs, all_haplotypes)

def get_subpopulation_haps(sample_IDs, all_haplotypes):

    sample_ID_stream = open(sample_IDs, 'r')
    pop_sample_IDs = []
    for line in sample_ID_stream:
        ID = line.strip().split()[0]
        pop_sample_IDs.append(ID)
    sample_ID_stream.close()


    all_hap_stream = open(all_haplotypes, 'r')
    for line in all_hap_stream:
        sample_hap = line.strip().split('\t')[0]
        root_ID = sample_hap.split('_')[0]
        if root_ID in pop_sample_IDs:
            print(sample_hap)
            #print(population_code, superpopulation_code)

    all_hap_stream.close()


if __name__ == '__main__':
    main()
