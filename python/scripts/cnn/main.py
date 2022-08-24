import basic_ds
import sys
import utils
sys.path.insert(1, '/home/sdp/precision-medicine/python/scripts/')
import basic_data_structures

import tensorflow as tf
import numpy as np

sample_IDs_file = sys.argv[1]
sample_encodings_file = sys.argv[2]
ID_distances_file = sys.argv[3]
ID_embeddings_file = sys.argv[4]
#ID_embeddings_file_2 = sys.argv[5]

# model parameters
num_epochs = 10
num_layers = 5
filter_size = 32


## rewriting my main_cnn.py script in-line like the example

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


## basic data structures
sample_IDs = basic_data_structures.get_sample_ID_list(sample_IDs_file)
sample_encodings = basic_data_structures.get_encoding_list(sample_encodings_file)


print('Building dataset...')
print('...sampling without replacement')
ds = basic_ds.sample_without_replacement(sample_IDs_file, sample_encodings_file,
                                         ID_distances_file,
                                         num_samples, num_variants)
print('Done.')
print()

print('Splitting into training and validation...')
train_dataset = ds.take(round(num_pairs * 0.8))
val_dataset = ds.skip(round(num_pairs * 0.8))

# getting size of the input encoding vectors
vector_size = num_variants

# building network
print('Building Base CNN...')
inputs = tf.keras.Input(shape=(num_variants, 1))

# using for loop to build layers
for l in range(num_layers):
    x = tf.keras.layers.Conv1D(filters=filter_size, kernel_size=3, activation="relu")(inputs)
    x = tf.keras.layers.MaxPool1D(pool_size=2, strides=2)(x)
    x = tf.keras.layers.BatchNormalization()(x)
    filter_size *= 2

x = tf.keras.layers.GlobalAveragePooling1D()(x)

outputs = tf.keras.layers.Dense(80)(x)

base_cnn = tf.keras.Model(inputs=inputs, outputs=outputs, name="base_cnn")
embedding = tf.keras.Model(base_cnn.input, outputs, name="Embedding")


# Distances layer computes the distance
# between two samples
class DistanceLayer(tf.keras.layers.Layer):
    """
    Computes distance between two vectors with simple euclidean distance
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def call(self, sample1, sample2):
        distance = tf.reduce_sum(tf.square(sample1 - sample2), -1)
        # distance = tf.norm(sample1 - sample2, ord='euclidean')
        return distance


sample1_input = tf.keras.Input(name="sample1", shape=vector_size)
sample2_input = tf.keras.Input(name="sample2", shape=vector_size)

predicted_distances = DistanceLayer()(
    embedding(sample1_input),
    embedding(sample2_input)
)
siamese_network = tf.keras.Model(inputs=[sample1_input, sample2_input],
                                 outputs=predicted_distances)


class SiameseModel(tf.keras.Model):
    """
    Computes loss using two embeddings proces by base CNN
    """
    def __init__(self, siamese_network):
        super(SiameseModel, self).__init__()
        self.siamese_network = siamese_network
        self.loss_tracker = tf.keras.metrics.Mean(name="loss")

    def call(self, data):
        sample1, sample2 = data[:2]
        # sample1 = tf.expand_dims(input=sample1, axis=-1)
        # sample2 = tf.expand_dims(input=sample2, axis=-1)

        return self.siamese_network(sample1, sample2)

    def train_step(self, data):
        """
        :param data: sample1, sample2
        :return:
        """
        sample_data = data[:2]
        target_data = data[2]

        with tf.GradientTape() as tape:
            loss = self._compute_loss(sample_data, target_data)

        gradients = tape.gradient(loss, self.siamese_network.trainable_weights)

        self.optimizer.apply_gradients(
            zip(gradients, self.siamese_network.trainable_weights)
        )

        self.loss_tracker.update_state(loss)
        return {"loss": self.loss_tracker.result()}

    def test_step(self, data):
        sample_data = data[:2]
        target_data = data[2]
        loss = self._compute_loss(sample_data, target_data)

        self.loss_tracker.update_state(loss)
        return {"loss": self.loss_tracker.result()}

    def _compute_loss(self, sample_data, target_distances):
        # predicted_distance = self.siamese_network(sample_data)
        predicted_distance = self.siamese_network(sample_data)

        loss = tf.keras.metrics.mean_squared_error(
            predicted_distance, target_distances
        )
        return loss

    @property
    def metrics(self):
        return [self.loss_tracker]


siamese_model = SiameseModel(siamese_network)

callback = tf.keras.callbacks.EarlyStopping(monitor='loss', min_delta=0)#mode="auto", patience=1)
siamese_model.compile(optimizer=tf.keras.optimizers.Adam(0.0001), run_eagerly=True, weighted_metrics=[])#, loss='loss')
siamese_model.fit(train_dataset, epochs=num_epochs, validation_data=val_dataset)

out_embeddings = open(ID_embeddings_file, 'w')

all_sample_embeddings = []
# iterate through all pairs in a batch
for batch in train_dataset:
    # get sample 1 and sample 2 in a pair
    sample1, sample2 = batch[:2]

    # embeddings
    sample1_embedding, sample2_embedding = (
        embedding(sample1),
        embedding(sample2)
    )
    
    for s in range(len(sample1_embedding.numpy())):
        all_sample_embeddings.append(sample1_embedding[s].numpy())
        all_sample_embeddings.append(sample2_embedding[s].numpy())

print("num embeddings: ", len(all_sample_embeddings))
for embedding_i in range(len(all_sample_embeddings)):
    curr_embedding = all_sample_embeddings[embedding_i]
    for f in curr_embedding:
        out_embeddings.write(str(f) + ' ')
    out_embeddings.write('\n')

out_embeddings.close()

    #for s in range(len(sample1_embedding)):
    #    np_embedding_1 = sample1_embedding[s].numpy()
    #    np_embedding_2 = sample2_embedding[s].numpy()
    #    for f in np_embedding_1:
    #        out_embeddings.write(str(f) + ' ')
    #    out_embeddings.write('\n')
    #    for f in np_embedding_2:
    #        out_embeddings.write(str(f) + ' ')
    #    out_embeddings.write('\n')

    #cosine_similarity = tf.metrics.CosineSimilarity()

    #encoding_similarity = cosine_similarity(sample1, sample2)
    #embedding_similarity = cosine_similarity(sample1_embedding, sample2_embedding)
    #print('Encoding CS: ', encoding_similarity.numpy(), 'Embedding CS: ', embedding_similarity.numpy(),
    #    'Delta: ', embedding_similarity.numpy() - encoding_similarity.numpy())


