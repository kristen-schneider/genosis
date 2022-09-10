import argparse
import read_ilash

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--ilash_matches')
    return parser.parse_args()

def main():
    args = get_args()
    ilash_pairs = read_ilash.make_ilash_pairs_dict(args.ilash_matches)

    all_pairs = []

    for sample1 in ilash_pairs:
        
        sample1_pairs = ilash_pairs[sample1]
        for sample2 in sample1_pairs:
            pair_data = sample1_pairs[sample2]
           
            #print(sample1, sample2, pair_data)
            all_pairs.append([sample1, sample2, pair_data[0], pair_data[1]])
    
            
            #if sample1 in sample2:
            #    print(sample1, sample2, pair_data)
    srt = sorted(all_pairs, key = lambda x: x[3], reverse=True)
    for s in srt: print(s)

    print(len(srt))
if __name__ == '__main__':
    main()
