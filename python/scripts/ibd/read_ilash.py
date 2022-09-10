from collections import defaultdict

def make_ilash_pairs_cumulative(ilash_file):
    '''
    reads output from ilash and generates a dictionary
    whos key is a sample name and whose value is another
    dictionary. the second dictionary has key as sample
    match (pair) and has value as tuple with (match_count, 
    cum_match_length).

    '''
    ilash_pairs_dict = defaultdict(dict)

    f = open(ilash_file)
    for line in f:
        L = line.strip().split()
        sample1_h = L[1].split('_')[0]
        sample2_h = L[3].split('_')[0]

        length_match = float(L[9])

        try:
            ilash_pairs_dict[sample1_h][sample2_h][0] += 1
            ilash_pairs_dict[sample1_h][sample2_h][1] += length_match
        except KeyError:
            try: 
                ilash_pairs_dict[sample1_h][sample2_h] = [1, length_match]
            except KeyError:
                ilash_pairs_dict[sample1_h] = {sample2_h: [1, length_match]}
    f.close()
    return ilash_pairs_dict


def make_ilash_pairs_individual(ilash_file):
    '''
    reads output from ilash and generates a dictionary
    whos key is a match index, and whose value is another
    dictionary. the second dictionary has key as sample
    pair [s1, s2] and has value as match length

    '''
    ilash_pairs_dict = defaultdict(dict)

    f = open(ilash_file)
    match_index = 0
    for line in f:
        L = line.strip().split()
        sample1_h = L[1]
        sample2_h = L[3]

        length_match = float(L[9])

        ilash_pairs_dict[match_index] = {sample1_h: {sample2_h: length_match}}
        match_index += 1
    f.close()
    return ilash_pairs_dict
