import argparse
import distance_calculations
import read_encoding

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--encoding_file')
    return parser.parse_args()

def main():
    args = get_args()
    encodings = read_encoding.read_encoding_file(args.encoding_file)

    query = encodings[0]
    print(query)
    database = encodings
    
    # hamming distance
    hamming_distances = []
    for v_i in database:
        v_hd = distance_calculations.hamming_distance(query, v_i)
        hamming_distances.append(v_hd)

    # normalized hamming distance
    norm_hamming_distances = []
    for hd in hamming_distances:
        n_hd = hd/len(query)
        norm_hamming_distances.append(n_hd)

    # euclidean distance
    euclidean_distances = []
    for v_i in database:
        v_ed = distance_calculations.euclidean_distance(query, v_i)
        euclidean_distances.append(v_ed)



if __name__ == '__main__':
    main()
