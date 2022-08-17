import tensorflow as tf
# tf.enable_eager_execution()

#def build_base_cnn(vector_size):
#    """
#    Builds a base CNN with provided vector_size
#    :param vector_size:
#    :return: base CNN (tf.keras.Model)
#    """
#    # Input Layer
#    inputs = tf.keras.Input(shape=(vector_size, 1))
#
#    # Convolution Layers
#    ## filters = number of filters
#    ## kernel_size = length of convolution window
#    x = tf.keras.layers.Conv1D(filters=6, kernel_size=3, activation="relu")(inputs)
#    x = tf.keras.layers.BatchNormalization()(x)
#    x = tf.keras.layers.Conv1D(filters=6, kernel_size=3, activation="relu")(x)
#    x = tf.keras.layers.BatchNormalization()(x)
#
#    x = tf.keras.layers.GlobalAveragePooling1D()(x)
#    # dense1 = tf.keras.layers.Dense(40, activation="relu", name="one")(inputs)
#    # dense2 = tf.keras.layers.Dense(20, activation="relu", name="two")(dense1)
#
#    # Outputs
#    outputs = tf.keras.layers.Dense(80)(x)
#    # outputs = tf.keras.layers.Dense(20)(x)
#
#    # base CNN model
#    base_cnn_model = tf.keras.Model(inputs=inputs, outputs=outputs, name="base_cnn")
#
#    return base_cnn_model

def build_embedding(vector_size):
    # Input Layer
    inputs = tf.keras.Input(shape=(vector_size, 1))

    # Convolution Layers
    ## filters = number of filters
    ## kernel_size = length of convolution window
    x = tf.keras.layers.Conv1D(filters=6, kernel_size=3, activation="relu")(inputs)
    x = tf.keras.layers.BatchNormalization()(x)
    x = tf.keras.layers.Conv1D(filters=6, kernel_size=3, activation="relu")(x)
    x = tf.keras.layers.BatchNormalization()(x)
    x = tf.keras.layers.Conv1D(filters=6, kernel_size=3, activation="relu")(x)
    x = tf.keras.layers.BatchNormalization()(x)

    x = tf.keras.layers.GlobalAveragePooling1D()(x)

    # Outputs
    outputs = tf.keras.layers.Dense(2500)(x)
    # outputs = tf.keras.layers.Dense(20)(x)

    # base CNN model
    base_cnn = tf.keras.Model(inputs=inputs, outputs=outputs, name="base_cnn")
    embedding = tf.keras.Model(base_cnn.input, outputs, name="Embedding")
    return embedding

def build_siamese_network(vector_size):
    """
    Builds siamese network
    :param vector_size: size of the sample inputs
    :return: siamese network with two inputs (s1, s2) and one output (predicted distance)
    """
    embedding = build_embedding(vector_size)

    sample1_input = tf.keras.Input(name="sample1", shape=vector_size)
    sample2_input = tf.keras.Input(name="sample2", shape=vector_size)

    predicted_distances = DistanceLayer()(
        embedding(sample1_input),
        embedding(sample2_input)
    )
    siamese_network = tf.keras.Model(inputs=[sample1_input, sample2_input],
                                     outputs=predicted_distances)
    return siamese_network

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

        # print(sample_data[0].numpy())
        # print(target_data.numpy())

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

