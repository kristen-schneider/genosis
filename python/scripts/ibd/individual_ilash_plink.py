import argparse
import matplotlib.pyplot as plt
import read_ilash
import read_plink

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--ilash_matches')
    parser.add_argument('--plink_dir')
    parser.add_argument('--out_file')
    return parser.parse_args()

def main():
    args = get_args()
    print('getting ilash data')
    ilash_pairs = read_ilash.make_ilash_pairs_individual(args.ilash_matches)

    index_ilash_pairs = []
    index_plink_pairs = dict()

    for index in ilash_pairs:

        print('getting plink dict for ilash match: ', index)
        plink_file = args.plink_dir + 'ibd.' + str(index) + '.genome'
        plink_pairs_dict = read_plink.make_plink_pairs_dict(plink_file)
        print('...plink file:', plink_file)
    
        for sample1 in ilash_pairs[index]:
            s1_plink_notation = sample1.split('_')[0]
            for sample2 in ilash_pairs[index][sample1]:
                s2_plink_notation = sample2.split('_')[0]
                pair_cm = sampl2 = ilash_pairs[index][sample1][sample2]
                #sample1 = pair[0]
                #sample2 = pair[1]
                #pair_cm = ilash_pairs[index][pair]

                index_ilash_pairs.append([sample1, sample2, pair_cm]) 
        
        
        try:
            index_plink_dist = plink_pairs_dict[s1_plink_notation][s2_plink_notation]
        except KeyError:
            try:
                index_plink_dist = plink_pairs_dict[s2_plink_notation][s1_plink_notation]
            except KeyError:
                print('iLASH reported same sample similarity' +
                '(could be different haplotypes).\n' +
                'Plink does not report same sample. Cannot compare:')
                print(s1_plink_notations, s2_plink_notation)
            

        try:
            index_plink_pairs[sample1][sample2].append(index_plink_dist)
        except KeyError: 
            try:
                index_plink_pairs[sample1] = {sample2: [index_plink_dist]}
            except KeyError:
                index_plink_pairs[sample1][sample2] = [index_plink_dist]


    ilash_sorted = sorted(all_ilash_pairs, key = lambda x: x[3], reverse=True)
    for s in ilash_sorted: print(s)

    write_ilash_plink(ilash_sorted, index_plink_pairs, args.out_file)

def write_ilash_plink(ilash_dict, plink_dict, out_file):

    plink_x = []
    ilash_y = []

    for data in range(len(ilash_dict)):
        sample1 = ilash_dict[data][0]
        sample2 = ilash_dict[data][1]
        cm_length = ilash_dict[data][3]

        ilash_y.append(cm_length)
        
        #try:
        plink_dist = plink_dict[sample1][sample2]
        #except KeyError:
        #    try:
        #        plink_dist = plink_dict[sample2][sample1]
        #    except KeyError:
        #        print('iLASH reported same sample similarity' + 
        #        '(could be different haplotypes).\n' +
        #        'Plink does not report same sample. Cannot compare:')
        #        print(sample1, sample2)

        plink_x.append(plink_dist)
        
    o = open(out_file, 'w')
    o.write('ilash plink\n')
    for i in range(len(ilash_y)):
        line = str(ilash_y[i]) + ' ' + str(plink_x[i]) + '\n'
        o.write(line)
    o.close()


if __name__ == '__main__':
    main()
