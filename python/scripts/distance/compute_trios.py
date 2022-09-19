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

    encoded_file = args.encoded_file
    encodings = read_encoding.read_encoding_file(encoded_file)
    query_hap = args.query+'_'+args.hap

    print("Query: ", query_hap)
    print("sample_ID dist")
    q_encoding = encodings[query_hap]

    gaps_allowed = 0
    for s in encodings:
        s_encoding = encodings[s]
        qs_kd = distance_calculations.kristen(q_encoding, s_encoding, gaps_allowed)
        print(s, qs_kd)


if __name__ == '__main__':
    main()
