from collections import defaultdict

def make_pibd_pairs_cumulative(pibd_file):
    '''
    reads output from phased IBD and generates a dictionary
    whos key is a sample name and whose value is another
    dictionary. the second dictionary has key as sample
    match (pair) and has value as tuple with (match_count, 
    cum_match_length).

    '''
    pibd_pairs_dict = defaultdict(dict)

    f = open(pibd_file)
    for line in f:
        L = line.strip().split()
        sample1_id = int(L[1])
        sample2_id = int(L[2])

        start_match = float(L[7])
        end_match = float(L[8])
        length_match = end_match - start_match

        try:
            pibd_pairs_dict[sample1_id][sample2_id][0] += 1
            pibd_pairs_dict[sample1_id][sample2_id][1] += length_match
        except KeyError:
            try: 
                pibd_pairs_dict[sample1_id][sample2_id] = [1, length_match]
            except KeyError:
                pibd_pairs_dict[sample1_id] = {sample2_id: [1, length_match]}
    f.close()
    return pibd_pairs_dict


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
