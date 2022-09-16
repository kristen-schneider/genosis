import argparse
import distance_calculations
import read_encoding

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--encoding_file')
    parser.add_argument('--query')
    parser.add_argument('--dad')
    parser.add_argument('--mom')
    parser.add_argument('--base')
    return parser.parse_args()

def main():
    args = get_args()

    print('Reading Encodings...')
    encodings = read_encoding.read_encoding_file(args.encoding_file)

    q = encodings[args.query]
    d = encodings[args.dad]
    m = encodings[args.mom]
    gaps_allowed = 0

    qq = distance_calculations.kristen(q, q, gaps_allowed)
    print(qq)
    qd = distance_calculations.kristen(q, d, gaps_allowed)
    print(qd)
    qm = distance_calculations.kristen(q, m, gaps_allowed)
    print(qm)

if __name__ == '__main__':
    main()
