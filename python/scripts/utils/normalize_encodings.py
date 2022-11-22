import argparse
import sys

sys.path.append('../distance/')
import distance_calculations
import read_encoding

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--encoded_file')
    return parser.parse_args()

def main():
    args = get_args()
    encoded_file = args.encoded_file
    encodings = read_encoding.read_encoding_file(encoded_file)
    for e in encodings:
        norm_encoding = normalize_vector(encodings[e])

        print(e, ' '.join([str(n) for n in norm_encoding]))

def normalize_vector(vector):
    '''
    sum all of the elements in the input vector, divide each element
    by the sum, and return the normalized vector.
    '''
    norm_vector = []
    v_sum = sum(vector)
    for v in vector:
        norm_vector.append(v/len(vector))
    return norm_vector

if __name__ == '__main__':
    main()
