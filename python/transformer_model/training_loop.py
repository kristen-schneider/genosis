# python libraries
import argparse
import numpy as np
import tensorflow as tf
    

# project scripts
#from attention import MultiHeadAttention
#from decoder import Decoder
#from encoder import Encoder
#from positional_encoding import PositionEmbeddingFixedWeights
from models import gt_transformer as gtt
from models import build_siamese_network
from models import SiameseModel
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

    num_variants = len(genotype_encodings_index[0])
    num_samples = training_IDs

    # convert input data to tensors
    gt_tensors = get_tensors(sample_names_index,
                    genotype_encodings_index,
                    training_IDs)
    train_ds = PairsDatasetUnsupervised(
        sample_id_filename=train_samples,
        genotype_filename=genotype_filename,
        shuffle=True,
        repeat=False,
        batch_size=batch_size,
    )

    # TODO:
    #   positional encoding with variable length input (basepair encoding)
    #   this does not work with variable length input
    pe_tensors = get_tensors(sample_names_index,
                    positional_encodings_index, 
                    training_IDs)

    # TODO: 
    #   sample pairs

    # Buildling base model (transformer positional encoding + encoder)
    # TODO:
    #   define input (pe input tensor)...?
    #   check base model for correctness...
    base_model = gtt(
            dim_val=512,
            num_encoder_layers=4,
            num_heads=8,
            dropout=0,
            activation='relu',
            layer_norm_epsilon=1e-05,
            kernel_initializer="glorot_uniform",
            bias_initializer="zeros",
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
