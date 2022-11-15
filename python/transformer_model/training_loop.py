# python libraries
import argparse
import numpy as np
import tensorflow as tf
    

# project scripts
#from attention import MultiHeadAttention
#from decoder import Decoder
#from encoder import Encoder
#from positional_encoding import PositionEmbeddingFixedWeights
from models import GTTransformer as gtt
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

    # TODO:
    #   positional encoding with variable length input (basepair encoding)
    #   
    #
    
    num_variants = len(genotype_encodings_index[0])

    # Buildling base model (transformer encoder)
    base_model = gtt(
                input_size = 64,
                out_seq_len = 58,
                dim_val = 512,
                num_encoder_layers = 4,
                num_heads = 8,
                dropout_encoder = 0.2,
                dropout_pos_enc = 0.1,
                dim_feedforward_encoder = 2048,
                activation = 'relu',
                )
    

    # Siamese network using base model
    siamese_net = build_siamese_network(
            base_model=base_model,
            vector_size=num_variants,
            )
    training_model = SiameseModel(siamese_net)
    training_model.compile(
            optimizer=tf.keras.optimizers.Adam(learning_rate=0.001),
            run_eargerly=True,
            )

    # TODO: pairs for training dataset
    training_model.fit(
            #training_IDs
            epochs=10,
            #...
            )



if __name__ == '__main__':
    main()
