import tensorflow as tf

class GTTransformer(tf.keras.Model):
    def __init__(self,
        input_size: int,
        out_seq_len: int=58,
        dim_val: int=512,
        num_encoder_layers: int=4,
        num_heads: int=8,
        dropout_encoder: float=0.2,
        dropout_pos_enc: float=0.1,
        dim_feedforward_encoder: int=2048,
        activation: str='relu',
        ):
        """
        """
        super().__init__() 

        # input layer        
        self.encoder_input_layer = tf.keras.layers.Dense(
            units=out_seq_length, 
            activation=activation 
            )
        # positional layer 
        # TODO: positional layer for variable-length input
        
        # attention layer (multihead)
        self.attention_layer = tf.keras.layers.MultiHeadAttention(
                num_heads,
                key_dim,)

        # encoder layer
        encoder_layer = tf.keras_nlp.layers.TransformerEncoder(
            dim_val,
            num_heads,
            dropout=dropout_encoder,
            activation=activation,
            layer_norm_epsilon=1e-05,
            kernel_initializer="glorot_uniform",
            bias_initializer="zeros",
            )
        
        # feed forward loop
        def feed_forward():
            x = self.encoder_input_layer(src)
            x = self.encoder_layer(x=x)
            return x



