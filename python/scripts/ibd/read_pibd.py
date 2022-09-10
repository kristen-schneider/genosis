from collections import defaultdict
import pysam

def make_pibd_pairs_cumulative(pibd_file, vcf_file):
    '''
    reads output from phased IBD and generates a dictionary
    whos key is a sample name and whose value is another
    dictionary. the second dictionary has key as sample
    match (pair) and has value as tuple with (match_count, 
    cum_match_length).

    '''
    sample_IDs = vcf_samples(vcf_file)
    pibd_pairs_dict = defaultdict(dict)

    f = open(pibd_file)
    header = f.readline()
    
    for line in f:
        L = line.strip().split()
        sample1_idx = int(L[2])
        sample1_name = sample_IDs[sample1_idx]
        sample2_idx = int(L[3])
        sample2_name = sample_IDs[sample2_idx]
    
        start_match = float(L[8])
        end_match = float(L[9])
        length_match = end_match - start_match

        try:
            pibd_pairs_dict[sample1_name][sample2_name][0] += 1
            pibd_pairs_dict[sample1_name][sample2_name][1] += length_match
        except KeyError:
            try: 
                pibd_pairs_dict[sample1_name][sample2_name] = [1, length_match]
            except KeyError:
                pibd_pairs_dict[sample1_name] = {sample2_name: [1, length_match]}
    f.close()
    return pibd_pairs_dict

def vcf_samples(vcf_file):
    '''
    reads vcf and retruens a list of sample IDs
    '''
    vcf_samples_list = []
    vcf_f = pysam.VariantFile(vcf_file)
    vcf_samples_list = list((vcf_f.header.samples))
    return vcf_samples_list






### NOT USING INDIVIDUAL ON PHASED IBD (for now) BECAUSE
### PHASED IBD REPORTS MANY IBD SEGMENTS
def make_pibd_pairs_individual(pibd_file):
    '''
    reads output from pibd and generates a dictionary
    whos key is a match index, and whose value is another
    dictionary. the second dictionary has key as sample1
    and has value as dictionary {sample2: match_length}

    '''
    pibd_pairs_dict = defaultdict(dict)

    f = open(pibd_file)
    match_index = 0
    for line in f:
        L = line.strip().split()
        sample1_h = L[1]
        sample2_h = L[3]

        length_match = float(L[9])

        pibd_pairs_dict[match_index] = {sample1_h: {sample2_h: length_match}}
        match_index += 1
    f.close()
    return pibd_pairs_dict
