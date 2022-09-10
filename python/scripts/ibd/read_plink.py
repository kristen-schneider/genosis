from collections import defaultdict

def make_plink_pairs_dict(plink_file, delim=' '):
    '''
    reads output from plink genome into a dictionary
    whos key is a sample name and whose value is another
    dictionary. the second dictionary has key as sample
    match (pair) and has value as the plink dist

    '''
    plink_pairs_dict = defaultdict(dict)

    f = open(plink_file, 'r')
    for line in f:
        L = line.strip().split()
        sample1 = L[1]
        sample2 = L[3]
        dist = L[11]

        plink_pairs_dict[sample1][sample2] = dist
    f.close()
    
    return plink_pairs_dict
