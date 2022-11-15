# python libraries
import argparse
import numpy as np
import tensorflow as tf
    

# project scripts
#from attention import MultiHeadAttention
#from decoder import Decoder
#from encoder import Encoder
#from positional_encoding import PositionEmbeddingFixedWeights
from models import GTTransformer
from utils import map_sample_names_to_index, \
                                map_genotype_encoding_to_index, \
                                map_positional_encoding_to_index, \
                                split_samples, \
                                get_tensors

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--sample_names')
    parser.add_argument('--sample_encodings')
    parser.add_argument('--sample_pos')
    return parser.parse_args()

    
def main():
    args = get_args()
    
    # loading sample names (sample_name: sample_index)
    sample_names_index = map_sample_names_to_index(args.sample_names)
    # loading genotype encodings (sample_index: genotype_encoding)
    genotype_encodings_index = map_genotype_encoding_to_index(args.sample_encodings)
    # loading positional encodings (sample_index: positional_encoding)
    positional_encodings_index = map_positional_encoding_to_index(args.sample_pos)
    # splitting IDs into training, testing, and validating sets
    training_IDs, testing_IDs, validating_IDs = split_samples(args.sample_names, .8)

    

if __name__ == '__main__':
    main()
