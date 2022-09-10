import argparse
import matplotlib.pyplot as plt
import read_pibd
import read_plink

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--pibd_matches')
    parser.add_argument('--plink_matches')
    parser.add_argument('--out_file')
    return parser.parse_args()

def main():
    args = get_args()
    print('getting pibd data')
    pibd_pairs = read_pibd.make_pibd_pairs_cumulative(args.pibd_matches)
    print('getting plink data')
    plink_pairs = read_plink.make_plink_pairs_dict(args.plink_matches)
    #print(plink_pairs)

    all_pibd_pairs = []
    for sample1 in pibd_pairs:
        sample1_pairs = pibd_pairs[sample1]
        for sample2 in sample1_pairs:
            pair_data = sample1_pairs[sample2] 
            all_pibd_pairs.append([sample1, sample2, pair_data[0], pair_data[1]]) 
    pibd_sorted = sorted(all_pibd_pairs, key = lambda x: x[3], reverse=True)
    for s in pibd_sorted: print(s)

    #write_pibd_plink(pibd_sorted, plink_pairs, args.out_file)

def write_pibd_plink(pibd_dict, plink_dict, out_file):

    plink_x = []
    pibd_y = []

    for data in pibd_dict:
        sample1 = data[0]
        sample2 = data[1]
        cm_length = data[3]

        pibd_y.append(cm_length)
        
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
    o.write('pibd plink\n')
    for i in range(len(pibd_y)):
        line = str(pibd_y[i]) + ' ' + str(plink_x[i]) + '\n'
        o.write(line)
    o.close()


if __name__ == '__main__':
    main()
