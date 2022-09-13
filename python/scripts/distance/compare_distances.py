import argparse
import distance_calculations
import read_encoding

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--encoding_file')
    parser.add_argument('--query')
    return parser.parse_args()

def main():
    args = get_args()
    encodings = read_encoding.read_encoding_file(args.encoding_file)

    query = encodings[args.query]
    database = encodings
    
    # hamming distance
    hamming_distances = dict()
    for v_i in database:
        v_hd = distance_calculations.hamming_distance(query, database[v_i])
        hamming_distances[v_i] = v_hd

    # normalized hamming distance
    norm_hamming_distances = dict()
    for hd in hamming_distances:
        n_hd = hamming_distances[hd]/len(query)
        norm_hamming_distances[hd] = n_hd

    # euclidean distance
    euclidean_distances = dict()
    for v_i in database:
        v_ed = distance_calculations.euclidean_distance(query, database[v_i])
        euclidean_distances[v_i] = v_ed



if __name__ == '__main__':
    main()
