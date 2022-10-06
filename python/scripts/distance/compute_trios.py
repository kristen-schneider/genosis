import argparse
import distance_calculations
import read_encoding

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--encoded_file')
    parser.add_argument('--hap')
    parser.add_argument('--query')
    return parser.parse_args()

def main():
    args = get_args()

    #print('Reading encodings...')
    encoded_file = args.encoded_file
    encodings = read_encoding.read_encoding_file(encoded_file)
    query_hap = args.query+'_'+args.hap

    print("Query: ", query_hap)
    print("sample_ID dist")
    q_encoding = encodings[query_hap]

    gaps_allowed = 0
    for s in encodings:
        s_encoding = encodings[s]
        s_base = '_'.join(s.split('_')[0:2])
        s_encoding_0 = encodings[s_base+'_'+str(0)]
        s_encoding_1 = encodings[s_base+'_'+str(1)]

        #qs_sv = distance_calculations.shared_variants(q_encoding, s_encoding)
        #qs_ed = distance_calculations.euclidean_distance(q_encoding, s_encoding)
        #qs_kd = distance_calculations.kristen(q_encoding, s_encoding, gaps_allowed)
        qs_rdp = distance_calculations.recombination_dp(q_encoding, s_encoding_0, s_encoding_1)
        print(s, qs_rdp)
        #if '_'+args.hap in s:
        #    s_encoding = encodings[s]
        #    #print(s_encoding)
        #    qs_ed = distance_calculations.euclidean_distance(q_encoding, s_encoding)
        #    #qs_kd = distance_calculations.kristen(q_encoding, s_encoding, gaps_allowed)
        #    print(s, qs_ed)


if __name__ == '__main__':
    main()
