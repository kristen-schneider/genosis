import tensorflow as tf

class GenotypeTransformer(tf.keras.Model):
    def __init__(self,
        input_size: int,
        out_seq_len: int=58,
        dim_val: int=512,
        n_encoder_layers: int=4,
        n_heads: int=8,
        dropout_encoder: float=0.2,
        dropout_pos_enc: float=0.1,
        dim_feedforward_encoder: int=2048,
        activation: str='relu',
        ):
        """
        """
        super().__init__() 


        #print("input_size is: {}".format(input_size))
        #print("dim_val is: {}".format(dim_val))

        # Creating the three linear layers needed for the model
        self.encoder_input_layer = tf.keras.layers.Dense(
            input_size, 
            activation=activation 
            )

        # Create positional encoder
        self.positional_encoding_layer = pe.PositionalEncoder(
            d_model=dim_val,
            dropout=dropout_pos_enc
            )
