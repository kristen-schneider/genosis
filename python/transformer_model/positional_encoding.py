"""
https://machinelearningmastery.com/a-gentle-introduction-to-positional-encoding-in-transformer-models-part-1
https://machinelearningmastery.com/the-transformer-positional-encoding-layer-in-keras-part-2/
"""

import numpy as np
import tensorflow as tf
from tensorflow.keras.layers import TextVectorization, Embedding, Layer

class PositionEmbeddingFixedWeights(Layer):
    def __init__(self, sequence_length, vocab_size, output_length, **kwargs):
        super(PositionEmbeddingFixedWeights, self).__init__(**kwargs)
        gt_embedding_matrix = self.get_position_encoding(vocab_size, output_length)
        position_embedding_matrix = self.get_position_encoding(sequence_length, output_length)
        self.gt_embedding_layer = Embedding(
            input_dim=vocab_size, output_dim=output_length,
            weights=[gt_embedding_matrix],
            trainable=False
        )
        self.position_embedding_layer = Embedding(
            input_dim=sequence_length, output_dim=output_length,
            weights=[position_embedding_matrix],
            trainable=False
        )

    def get_position_encoding(self, seq_len, d, n=10000):
        P = np.zeros((seq_len, d))
        for k in range(seq_len):
            for i in np.arange(int(d / 2)):
                denominator = np.power(n, 2 * i / d)
                P[k, 2 * i] = np.sin(k / denominator)
                P[k, 2 * i + 1] = np.cos(k / denominator)
        return P

    def call(self, inputs):
        position_indices = tf.range(tf.shape(inputs)[-1])
        embedded_gts = self.gt_embedding_layer(inputs)
        embedded_indices = self.position_embedding_layer(position_indices)
        return embedded_gts + embedded_indices
