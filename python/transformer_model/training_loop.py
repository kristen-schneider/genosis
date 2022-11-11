# python libraries
import argparse
import numpy as np
import tensorflow as tf

# project scripts
from attention import MultiHeadAttention
from encoder import Encoder
from positional_encoding import PositionEmbeddingFixedWeights
from utils import map_sample_names_to_index, \
                                map_genotype_encoding_to_index, \
                                split_samples, \
                                get_genotype_tensors

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--sample_names')
    parser.add_argument('--sample_encodings')
    return parser.parse_args()

def main():
    args = get_args()

    # loading sample names (sample_name: sample_index)
    sample_names_index = map_sample_names_to_index(args.sample_names)
    # loading genotype encodings (sample_index: genotype_encoding)
    genotype_encodings_index = map_genotype_encoding_to_index(args.sample_encodings)
    # splitting IDs into training, testing, and validating sets
    training_IDs, testing_IDs, validating_IDs = split_samples(args.sample_names, .8)

    # TODO: check these variables for correctness
    sequence_length = len(genotype_encodings_index[0])
    vocab_size = 4
    output_length = sequence_length
    # converting genotypes to tensors
    gt_tensors = get_genotype_tensors(sample_names_index,
                       genotype_encodings_index,
                       training_IDs)
    # converting tensors to fixed weight embeddings for input to positional encoding
    fixed_weights_embedding_layer = PositionEmbeddingFixedWeights(sequence_length,
                                                                  vocab_size,
                                                                  output_length)
    fixed_embedding = fixed_weights_embedding_layer(gt_tensors)

    # TODO: check these variables for correctness
    h = 8  # Number of self-attention heads
    d_k = 64  # Dimensionality of the linearly projected queries and keys
    d_v = 64  # Dimensionality of the linearly projected values
    d_model = 512  # Dimensionality of the model sub-layers' outputs
    batch_size = 64  # Batch size from the training process
    # TODO: change these to match data from genotype tensors
    queries = np.random.random((batch_size, sequence_length, d_k))
    keys = np.random.random((batch_size, sequence_length, d_k))
    values = np.random.random((batch_size, sequence_length, d_v))

    # mulit-head attention
    multihead_attention = MultiHeadAttention(h, d_k, d_v, d_model)
    print(multihead_attention(queries, keys, values))
    # TODO: check these variables for correctness
    n = 6  # Number of layers in the encoder stack
    d_ff = 2048  # Dimensionality of the inner fully connected layer
    dropout_rate = 0.1  # Frequency of dropping the input units in the dropout layers

    input_seq = np.random.random((batch_size, sequence_length))

    encoder = Encoder(vocab_size, sequence_length, h, d_k, d_v, d_model, d_ff, n, dropout_rate)
    print(encoder(input_seq, None, True))


if __name__ == '__main__':
    main()


