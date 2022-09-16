import argparse
import distance_calculations
import read_encoding
import sys

sys.path.insert(0, '../')
from ibd import read_plink

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--encoding_file')
    parser.add_argument('--query')
    parser.add_argument('--plink')
    parser.add_argument('--distances_dir')
    parser.add_argument('--base')
    return parser.parse_args()

def main():
    args = get_args()

    print('Reading Plink file...')
    plink_dict = read_plink.make_plink_pairs_dict(args.plink)               
    plink_query = args.query.split('_')[0]
    o = open(args.distances_dir + args.base + '.plink', 'w')
    o.write('query: ' + plink_query + '\n')
    o.write('sample_ID plink_distance\n')
    for sample_i in plink_dict[plink_query]:
        o.write(sample_i + ' ' + plink_dict[plink_query][sample_i] + '\n')
    o.close()

    print('Reading Encodings...')
    encodings = read_encoding.read_encoding_file(args.encoding_file)

    query = encodings[args.query]
    database = encodings
    
    print('...computing hamming distance')
    # hamming distance
    o = open(args.distances_dir + args.base + '.hammingD', 'w')
    hamming_distances = dict()
    for v_i in database:
        v_hd = distance_calculations.hamming_distance(query, database[v_i])
        hamming_distances[v_i] = v_hd
    o.write('query: ' + args.query + '\n')
    o.write('sample_ID hamming_distance\n')
    for d in hamming_distances:
        o.write(d + ' ' + str(hamming_distances[d])+ '\n')
    o.close()

    print('...computing hamming distance, normalized')
    # normalized hamming distance
    o = open(args.distances_dir + args.base + '.hammingDN', 'w')
    norm_hamming_distances = dict()
    for hd in hamming_distances:
        n_hd = hamming_distances[hd]/len(query)
        norm_hamming_distances[hd] = n_hd
    o.write('query: ' + args.query + '\n')
    o.write('sample_ID hamming_distance_normalized\n')
    for d in norm_hamming_distances:
        o.write(d + ' ' + str(norm_hamming_distances[d])+ '\n')
    o.close()

    print('...computing euclidean distance')
    # euclidean distance
    o = open(args.distances_dir + args.base + '.euclideanD', 'w')
    euclidean_distances = dict()
    for v_i in database:
        v_ed = distance_calculations.euclidean_distance(query, database[v_i])
        euclidean_distances[v_i] = v_ed
    o.write('query: ' + args.query + '\n')
    o.write('sample_ID euclidean_distance\n')
    for d in euclidean_distances:
        o.write(d + ' ' + str(euclidean_distances[d])+ '\n')
    o.close()

    print('...computing new distance')
    # new distance
    gap_allowed = 0
    o = open(args.distances_dir + args.base + '.newD', 'w')
    new_distances = dict()
    for v_i in database:
        v_nd = distance_calculations.kristen(query, database[v_i], gap_allowed)
        new_distances[v_i] = v_nd
    o.write('query: ' + args.query + '\n')
    o.write('sample_ID new_distance\n')
    for d in new_distances:
        o.write(d + ' ' + str(new_distances[d])+ '\n')
    o.close()


if __name__ == '__main__':
    main()
