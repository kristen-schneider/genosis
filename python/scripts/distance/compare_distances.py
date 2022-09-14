import argparse
import distance_calculations
import read_encoding

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--encoding_file')
    parser.add_argument('--query')
    parser.add_argument('--distances_dir')
    return parser.parse_args()

def main():
    args = get_args()
    print('Reading Encodings...')
    encodings = read_encoding.read_encoding_file(args.encoding_file)

    query = encodings[args.query]
    database = encodings
    
    print('...computing hamming distance')
    # hamming distance
    o = open(args.distances_dir + 'chr8.seg.0.hammingD', 'w')
    hamming_distances = dict()
    for v_i in database:
        v_hd = distance_calculations.hamming_distance(query, database[v_i])
        hamming_distances[v_i] = v_hd
    o.write('query: ' + args.query + '\n')
    o.write('sample hamming_distance\n')
    for d in hamming_distances:
        o.write(d + ' ' + str(hamming_distances[d])+ '\n')
    o.close()

    print('...computing hamming distance, normalized')
    # normalized hamming distance
    o = open(args.distances_dir + 'chr8.seg.0.hammingDN', 'w')
    norm_hamming_distances = dict()
    for hd in hamming_distances:
        n_hd = hamming_distances[hd]/len(query)
        norm_hamming_distances[hd] = n_hd
    o.write('query: ' + args.query + '\n')
    o.write('sample hamming_distance_normalized\n')
    for d in norm_hamming_distances:
        o.write(d + ' ' + str(norm_hamming_distances[d])+ '\n')
    o.close()

    print('...computing euclidean distance')
    # euclidean distance
    o = open(args.distances_dir + 'chr8.seg.0.euclideanD', 'w')
    euclidean_distances = dict()
    for v_i in database:
        v_ed = distance_calculations.euclidean_distance(query, database[v_i])
        euclidean_distances[v_i] = v_ed
    o.write('query: ' + args.query + '\n')
    o.write('sample euclidean_distance\n')
    for d in euclidean_distances:
        o.write(d + ' ' + str(euclidean_distances[d])+ '\n')
    o.close()



if __name__ == '__main__':
    main()
