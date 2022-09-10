import argparse
import matplotlib.pyplot as plt
import read_ilash
import read_plink

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--ilash_matches')
    parser.add_argument('--plink_matches')
    parser.add_argument('--out_file')
    return parser.parse_args()

def main():
    args = get_args()
    print('getting ilash data')
    ilash_pairs = read_ilash.make_ilash_pairs_cumulative(args.ilash_matches)
    print('getting plink data')
    plink_pairs = read_plink.make_plink_pairs_dict(args.plink_matches)
    #print(plink_pairs)

    all_ilash_pairs = []
    for sample1 in ilash_pairs:
        sample1_pairs = ilash_pairs[sample1]
        for sample2 in sample1_pairs:
            pair_data = sample1_pairs[sample2] 
            all_ilash_pairs.append([sample1, sample2, pair_data[0], pair_data[1]]) 
    ilash_sorted = sorted(all_ilash_pairs, key = lambda x: x[3], reverse=True)
    for s in ilash_sorted: print(s)

    write_ilash_plink(ilash_sorted, plink_pairs, args.out_file)

def write_ilash_plink(ilash_dict, plink_dict, out_file):

    plink_x = []
    ilash_y = []

    for data in ilash_dict:
        sample1 = data[0]
        sample2 = data[1]
        cm_length = data[3]

        ilash_y.append(cm_length)
        
        try:
            plink_dist = plink_dict[sample1][sample2]
        except KeyError:
            try:
                plink_dist = plink_dict[sample2][sample1]
            except KeyError:
                print('iLASH reported same sample similarity' + 
                '(could be different haplotypes).\n' +
                'Plink does not report same sample. Cannot compare:')
                print(sample1, sample2)

        plink_x.append(plink_dist)
        
    o = open(out_file, 'w')
    o.write('ilash plink\n')
    for i in range(len(ilash_y)):
        line = str(ilash_y[i]) + ' ' + str(plink_x[i]) + '\n'
        o.write(line)
    o.close()


if __name__ == '__main__':
    main()
