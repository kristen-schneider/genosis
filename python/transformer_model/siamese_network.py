import tensorflow as tf


# def conv1d_block(
#     inputs: tf.Tensor,
#     filters: int,
#     kernel_size: int,
#     strides: int = 1,
#     activation: str = "relu",
# ) -> tf.Tensor:
#     """
#     Creates a conv1d block.
#     :param inputs: The input tensor.
#     :param filters: The number of filters.
#     :param kernel_size: The kernel size.
#     :param strides: stride length.
#     :return: The output tensor.
#     """
#     # x = tf.keras.layers.Conv1D(filters, kernel_size, strides=strides, padding="same")(
#     #     inputs
#     # )
#     x = tf.keras.layers.Conv1D(
#         filters,
#         kernel_size,
#         strides=strides,
#         padding="same",
#         kernel_regularizer=tf.keras.regularizers.L2(0.001),
#     )(inputs)
#     x = tf.keras.layers.BatchNormalization()(x)
#     x = tf.keras.layers.Activation(activation)(x)
#     return x
#
#
# def residual1d_block(
#     inputs: tf.Tensor,
#     filters: int,
#     kernel_size: int,
#     strides: int = 1,
#     activation: str = "relu",
# ) -> tf.Tensor:
#     """
#     Creates a residual conv1d block.
#     :param inputs: The input tensor.
#     :param filters: The number of filters.
#     :param kernel_size: The kernel size.
#     :param strides: stride length.
#     :return: The output tensor.
#     """
#
#     x = tf.keras.layers.Conv1D(
#         filters,
#         kernel_size,
#         strides=strides,
#         padding="same",
#         kernel_regularizer=tf.keras.regularizers.L2(0.001),
#     )(inputs)
#     x = tf.keras.layers.BatchNormalization()(x)
#     x = tf.keras.layers.Activation(activation)(x)
#     x = tf.keras.layers.Conv1D(
#         filters,
#         kernel_size,
#         strides=strides,
#         padding="same",
#         kernel_regularizer=tf.keras.regularizers.L2(0.001),
#     )(x)
#     x = tf.keras.layers.BatchNormalization()(x)
#     x = tf.keras.layers.Add()([x, inputs])
#     x = tf.keras.layers.Activation(activation)(x)
#     return x
#
#
# # define a resnet model using the above functions
# def resnet_model(
#     shape: tuple[int, int],  # length, channels
#     # inputs: tf.Tensor,
#     # filters: int,
#     # kernel_size: int,
#     # strides: int = 1,
#     activation: str = "relu",
#     n_residual_blocks: int = 3,
#     n_dense_blocks: int = 1,
#     dense_size: int = 64,
# ) -> tf.keras.Model:
#     """
#     Creates a residual conv1d block.
#     :param inputs: The input tensor.
#     :param filters: The number of filters.
#     :param kernel_size: The kernel size.
#     :param strides: stride length.
#     :return: The output tensor.
#     """
#     inputs = tf.keras.Input(shape=shape)
#     x = conv1d_block(inputs, filters=32, kernel_size=7)
#     x = tf.keras.layers.MaxPooling1D(pool_size=3, strides=1)(x)
#
#     for _ in range(n_residual_blocks):
#         x = residual1d_block(x, filters=32, kernel_size=3, activation=activation)
#
#     x = tf.keras.layers.GlobalAveragePooling1D()(x)
#     for _ in range(n_dense_blocks):
#         x = tf.keras.layers.Dense(dense_size, activation=activation)(x)
#     return tf.keras.Model(inputs=inputs, outputs=x)


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


def build_siamese_network(base_model: tf.keras.Model, vector_size: int):
    """
    Builds siamese network
    :param vector_size: size of the sample inputs
    :return: siamese network with two inputs (s1, s2) and one output (predicted distance)
    """

    # TODO allow for channels dimension
    sample1_input = tf.keras.Input(name="sample1", shape=vector_size)
    sample2_input = tf.keras.Input(name="sample2", shape=vector_size)

    # predicted_distances = DistanceLayer()(
    predicted_distances = tf.keras.layers.Dot(axes=-1, normalize=True)(
        [base_model(sample1_input), base_model(sample2_input)]
    )
    siamese_network = tf.keras.Model(
        inputs=[sample1_input, sample2_input], outputs=predicted_distances
    )
    return siamese_network


class SiameseModel(tf.keras.Model):
    """
    Wraps siamese network for training
    """

    def __init__(self, base_model: tf.keras.Model, vector_size: int):
        super(SiameseModel, self).__init__()
        self.base_model = base_model  # keep reference to base model so we can save it
        self.siamese_network = build_siamese_network(self.base_model, vector_size)
        self.loss_tracker = tf.keras.metrics.Mean(name="loss")

    def call(self, data):
        sample1, sample2 = data[:2]
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
        predicted_distance = self.siamese_network(sample_data)
        # loss = tf.keras.metrics.mean_squared_error(predicted_distance, target_distances)
        loss = tf.keras.losses.Huber()(predicted_distance, target_distances)
        return loss

    @property
    def metrics(self):
        return [self.loss_tracker]

    def save(self, *args, **kwargs):
        """
        Custom save method to only save the base model.
        """
        self.base_model.save(*args, **kwargs)


class SiameseNN(tf.keras.Model):
    """
    Build siamese network with simSiam style architecture
    """

    def __init__(
        self,
        input_dim: tuple[int, int],  # length, channels
        embedding_dim: int,
        activation: str = "relu",
        n_residual_blocks: int = 3,  # number of res blocks in base model
        n_dense_blocks: int = 1,  # number of dense layers in base model
    ):
        super(SiameseNN, self).__init__()

        self.input_dim = input_dim
        self.embedding_dim = embedding_dim
        self.activation = activation

        # TODO make this a transformer model
        self.base_model = resnet_model(
            shape=input_dim,
            activation=activation,
            n_residual_blocks=n_residual_blocks,
            n_dense_blocks=n_dense_blocks,
        )
        self.siamese_network = self._build_siamese_network()
        self.loss_tracker = tf.keras.metrics.Mean(name="loss")

    # TODO implement method to build siamese network
    def _build_siamese_network(self):
        """
        Builds siamese network
        :param vector_size: size of the sample inputs
        :return: siamese network with two inputs (s1, s2) and one output (predicted distance)
        """

        # TODO allow general input shape tuple
        in1 = tf.keras.Input(name="sample1", shape=self.input_dim[0])
        in2 = tf.keras.Input(name="sample2", shape=self.input_dim[0])

        # Compute embeddings. Only track the gradient through 1 leg of the network
        x1 = self.base_model(in1)
        # x2 = tf.stop_gradient(self.base_model(in2))
        x2 = self.base_model(in2)

        pred_sim = tf.keras.layers.Dot(axes=-1, normalize=True)([x1, x2])
        siamese_network = tf.keras.Model(inputs=[in1, in2], outputs=pred_sim)
        return siamese_network

    def call(self, data):
        return self.siamese_network(data)

    def train_step(self, data):
        with tf.GradientTape() as tape:
            loss = self._compute_loss(data)
        gradients = tape.gradient(loss, self.siamese_network.trainable_weights)
        self.optimizer.apply_gradients(
            zip(gradients, self.siamese_network.trainable_weights)
        )
        self.loss_tracker.update_state(loss)
        return {"loss": self.loss_tracker.result()}


    def test_step(self, data):
        loss = self._compute_loss(data)
        self.loss_tracker.update_state(loss)
        return {"loss": self.loss_tracker.result()}

    def _compute_loss(self, data):
        pred_sim = self.call(data)

        # we are using the full dimension similarity as
        # the ground truth, and computing it on the fly
        true_sim = tf.stop_gradient(
            tf.keras.layers.Dot(axes=-1, normalize=True)([data[0], data[1]])
        )
        return tf.keras.losses.MSE(true_sim, pred_sim)

    @property
    def metrics(self):
        return [self.loss_tracker]

    def save(self, *args, **kwargs):
        """
        Custom save method to only save the base model.
        """
        self.base_model.save(*args, **kwargs)
