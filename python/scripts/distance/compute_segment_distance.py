import argparse
import distance_calculations
import read_encoding

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--encoded_file')
    parser.add_argument('--query_file')
    return parser.parse_args()

def main():
    args = get_args()
    encoded_file = args.encoded_file
    encodings = read_encoding.read_encoding_file(encoded_file)
    
    queries = read_queries(args.query_file)
    haps = [0,1]
    for q in queries:
        for hap in haps:
            query_hap = q+str('_')+str(hap)
            print("Query: ", query_hap)
            q_encoding = encodings[query_hap]
            for s in encodings:
                s_encoding = encodings[s]
                #dist = distance_calculations.euclidean_distance(q_encoding, s_encoding)
                dist = distance_calculations.cosine_similarity(q_encoding, s_encoding)
                print(s, dist)
            print()

def read_queries(query_file):
    '''
    read a file of sample IDs for queries
    and return a list
    '''
    queries = []
    f = open(query_file, 'r')
    for line in f:
        q = line.strip()
        queries.append(q)
    f.close()
    return queries

if __name__ == '__main__':
    main()
