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

    print('-Loading sample names...')
    sample_names_index = map_sample_names_to_index(args.sample_names)
    print('-Loading genotype encodings...')
    genotype_encodings_index = map_genotype_encoding_to_index(args.sample_encodings)
    print('-Splitting samples into training, testing, and validating set...')
    training_IDs, testing_IDs, validating_IDs = split_samples(args.sample_names, .8)

    sequence_length = len(genotype_encodings_index[0])
    vocab_size = 4
    output_length = sequence_length
    print('-Converting genotypes to tensors...')
    gt_tensors = get_genotype_tensors(sample_names_index,
                       genotype_encodings_index,
                       training_IDs)
    print('-Converting tensors to fixed weight embeddings...')
    fixed_weights_embedding_layer = PositionEmbeddingFixedWeights(sequence_length,
                                                                  vocab_size,
                                                                  output_length)
    fixed_embedding = fixed_weights_embedding_layer(gt_tensors)
    print('-Testing encoder...')

    h = 8  # Number of self-attention heads
    d_k = 64  # Dimensionality of the linearly projected queries and keys
    d_v = 64  # Dimensionality of the linearly projected values
    d_model = 512  # Dimensionality of the model sub-layers' outputs
    batch_size = 64  # Batch size from the training process

    queries = np.random.random((batch_size, sequence_length, d_k))
    keys = np.random.random((batch_size, sequence_length, d_k))
    values = np.random.random((batch_size, sequence_length, d_v))

    multihead_attention = MultiHeadAttention(h, d_k, d_v, d_model)
    print(multihead_attention(queries, keys, values))

    n = 6  # Number of layers in the encoder stack
    d_ff = 2048  # Dimensionality of the inner fully connected layer
    dropout_rate = 0.1  # Frequency of dropping the input units in the dropout layers

    input_seq = np.random.random((batch_size, sequence_length))

    encoder = Encoder(vocab_size, sequence_length, h, d_k, d_v, d_model, d_ff, n, dropout_rate)
    print(encoder(input_seq, None, True))


if __name__ == '__main__':
    main()


