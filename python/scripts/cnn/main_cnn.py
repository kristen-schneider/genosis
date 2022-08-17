"""
Title: Similarity between two vectors for IBD
Authors: Kristen Schneider
Date created: 2022/07/20
Last modified: 2022/07/20
Description: Inspired by [code](https://keras.io/examples/vision/siamese_network/) written by [Hazem Essam](https://twitter.com/hazemessamm) and [Santiago L. Valdarrama](https://twitter.com/svpino)
"""

import basic_ds
import model
import sys
import utils

import tensorflow as tf

sample_IDs_file = sys.argv[1]
sample_encodings_file = sys.argv[2]
ID_distances_file = sys.argv[3]
#encoding_distances_file = sys.argv[4]
# ds_out_dir = sys.argv[4]
#tf_records_dir = sys.argv[4]

def main():

    # find dimensions of incoming data
    #   number of variants = size of input vector
    #   number of distances = number of pairwise distances computed
    [num_variants, num_samples] = utils.get_dimensions(sample_encodings_file)
    num_pairs = utils.count_pairs(num_samples)
    print('Opened file: ', sample_encodings_file)
    print('...Number of Variants: ', num_variants)
    print('...Number of Samples: ', num_samples)
    print('...Number of Pairs: ', num_pairs)
    print()
    
    print('Building dataset...')
    # ds = basic_ds.build_dataset_from_file(CNN_input_file)
    ds = basic_ds.sample_without_replacement(sample_IDs_file, sample_encodings_file,
                                             ID_distances_file,
                                             num_samples, num_variants)
    # ds = basic_ds.build_dataset_from_file(encoding_distances_file)
    # print('Writing dataset to zip file...')
    # tf.data.experimental.save(
    #     ds, ds_output_file, compression='GZIP'
    # )
    # load dataset from file
    # ds = tf.data.experimental.load(
    #     ds_output_file, element_spec=None, compression='GZIP', reader_func=None
    # )
    print('Done.')
    print()

    print('Splitting into training and validation...')
    train_dataset = ds.take(round(num_pairs * 0.8))
    val_dataset = ds.skip(round(num_pairs * 0.8))

    # val_dataset = val_dataset.batch(32, drop_remainder=False)
    # val_dataset = val_dataset.prefetch(tf.data.AUTOTUNE)

    # getting size of the input encoding vectors
    vector_size = num_variants

    # building network
    print('Building Network...')
    siamese_network = model.build_siamese_network(vector_size)

    # building model
    print('Building Model...')
    siamese_model = model.SiameseModel(siamese_network)
    siamese_model.compile(optimizer=tf.keras.optimizers.Adam(0.0001), run_eagerly=True)

    # training model
    print('Training Model...')
    siamese_model.fit(train_dataset, epochs=1, validation_data=val_dataset)

    # getting size of the input encoding vectors
    vector_size = num_variants

    print('Embeddings...')
    embedding = model.build_embedding(vector_size)
    
    # getting size of the input encoding vectors
    vector_size = num_variants

    print('Embeddings...')
    embedding = model.build_embedding(vector_size)
    
    count = 0
    for batch in train_dataset:
        count += 1
        sample1, sample2 = batch[:2]
        sample1_embedding, sample2_embedding = (
            embedding(sample1),
            embedding(sample2)
        )
        print(sample1, sample1_embedding)
        print(sample2, sample2_embedding)
    print(count)

    #sample = next(iter(train_dataset))
    #
    #sample1, sample2 = sample[:2]
    #sample1_embedding, sample2_embedding = (
    #    embedding(sample1),
    #    embedding(sample2)
    #)
    #
    #cosine_similarity = tf.metrics.CosineSimilarity()
    #similarity = cosine_similarity(sample1_embedding, sample2_embedding)
    #print("Similarity:", similarity.numpy())

if __name__ == '__main__':
    main()
