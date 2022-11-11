# python libraries
import argparse
from keras.losses import sparse_categorical_crossentropy
from model import TransformerModel
import numpy as np
from prepare_dataset import PrepareDataset
import tensorflow as tf
from tensorflow.keras.optimizers import Adam
from tensorflow.keras.optimizers.schedules import LearningRateSchedule
from tensorflow.keras.metrics import Mean
from tensorflow import data, train, math, reduce_sum, cast, equal, argmax, float32, GradientTape, TensorSpec, function, int64
from time import time
    

# project scripts
from attention import MultiHeadAttention
from decoder import Decoder
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

    # loading sample names (sample_name: sample_index)
    sample_names_index = map_sample_names_to_index(args.sample_names)
    # loading genotype encodings (sample_index: genotype_encoding)
    genotype_encodings_index = map_genotype_encoding_to_index(args.sample_encodings)
    # splitting IDs into training, testing, and validating sets
    training_IDs, testing_IDs, validating_IDs = split_samples(args.sample_names, .8)

    # TODO: check these variables for correctness
    sequence_length = len(genotype_encodings_index[0])
    vocab_size = 4
    output_length = sequence_length
    # converting genotypes to tensors
    gt_tensors = get_genotype_tensors(sample_names_index,
                       genotype_encodings_index,
                       training_IDs)
    # converting tensors to fixed weight embeddings for input to positional encoding
    fixed_weights_embedding_layer = PositionEmbeddingFixedWeights(sequence_length,
                                                                  vocab_size,
                                                                  output_length)
    fixed_embedding = fixed_weights_embedding_layer(gt_tensors)

    # TODO: check these variables for correctness
    h = 8  # Number of self-attention heads
    d_k = 64  # Dimensionality of the linearly projected queries and keys
    d_v = 64  # Dimensionality of the linearly projected values
    d_model = 512  # Dimensionality of the model sub-layers' outputs
    batch_size = 64  # Batch size from the training process
    # TODO: change these to match data from genotype tensors
    queries = np.random.random((batch_size, sequence_length, d_k))
    keys = np.random.random((batch_size, sequence_length, d_k))
    values = np.random.random((batch_size, sequence_length, d_v))

    # mulit-head attention
    multihead_attention = MultiHeadAttention(h, d_k, d_v, d_model)
    print(multihead_attention(queries, keys, values))
    # TODO: check these variables for correctness
    n = 6  # Number of layers in the encoder stack
    d_ff = 2048  # Dimensionality of the inner fully connected layer
    dropout_rate = 0.1  # Frequency of dropping the input units in the dropout layers

    input_seq = np.random.random((batch_size, sequence_length))

    encoder = Encoder(vocab_size, sequence_length, h, d_k, d_v, d_model, d_ff, n, dropout_rate)
    print(encoder(input_seq, None, True))


def procedures():
    """
    https://machinelearningmastery.com/training-the-transformer-model
    """
    # Define the model parameters
    h = 8  # Number of self-attention heads
    d_k = 64  # Dimensionality of the linearly projected queries and keys
    d_v = 64  # Dimensionality of the linearly projected values
    d_model = 512  # Dimensionality of model layers' outputs
    d_ff = 2048  # Dimensionality of the inner fully connected layer
    n = 6  # Number of layers in the encoder stack
    
    # Define the training parameters
    epochs = 2
    batch_size = 64
    beta_1 = 0.9
    beta_2 = 0.98
    epsilon = 1e-9
    dropout_rate = 0.1
    
    
    # Implementing a learning rate scheduler
    class LRScheduler(LearningRateSchedule):
        def __init__(self, d_model, warmup_steps=4000, **kwargs):
            super(LRScheduler, self).__init__(**kwargs)
    
            self.d_model = cast(d_model, float32)
            self.warmup_steps = warmup_steps
    
        def __call__(self, step_num):
    
            # Linearly increasing the learning rate for the first warmup_steps, and decreasing it thereafter
            arg1 = step_num ** -0.5
            arg2 = step_num * (self.warmup_steps ** -1.5)
    
            return (self.d_model ** -0.5) * math.minimum(arg1, arg2)
    
    
    # Instantiate an Adam optimizer
    optimizer = Adam(LRScheduler(d_model), beta_1, beta_2, epsilon)
   
    # TODO: CHANGE THIS FOR OUR GENOTYPE ENCODING INPUT
    # Prepare the training and test splits of the dataset
    dataset = PrepareDataset()
    trainX, trainY, train_orig, enc_seq_length, dec_seq_length, enc_vocab_size, dec_vocab_size = dataset('english-german-both.pkl')
    
    # Prepare the dataset batches
    train_dataset = data.Dataset.from_tensor_slices((trainX, trainY))
    train_dataset = train_dataset.batch(batch_size)
    
    # Create model
    training_model = TransformerModel(enc_vocab_size, dec_vocab_size, enc_seq_length, dec_seq_length, h, d_k, d_v, d_model, d_ff, n, dropout_rate)
    
    
    # Defining the loss function
    def loss_fcn(target, prediction):
        # Create mask so that the zero padding values are not included in the computation of loss
        padding_mask = math.logical_not(equal(target, 0))
        padding_mask = cast(padding_mask, float32)
    
        # Compute a sparse categorical cross-entropy loss on the unmasked values
        loss = sparse_categorical_crossentropy(target, prediction, from_logits=True) * padding_mask
    
        # Compute the mean loss over the unmasked values
        return reduce_sum(loss) / reduce_sum(padding_mask)
    
    
    # Defining the accuracy function
    def accuracy_fcn(target, prediction):
        # Create mask so that the zero padding values are not included in the computation of accuracy
        padding_mask = math.logical_not(equal(target, 0))
    
        # Find equal prediction and target values, and apply the padding mask
        accuracy = equal(target, argmax(prediction, axis=2))
        accuracy = math.logical_and(padding_mask, accuracy)
    
        # Cast the True/False values to 32-bit-precision floating-point numbers
        padding_mask = cast(padding_mask, float32)
        accuracy = cast(accuracy, float32)
    
        # Compute the mean accuracy over the unmasked values
        return reduce_sum(accuracy) / reduce_sum(padding_mask)
    
    
    # Include metrics monitoring
    train_loss = Mean(name='train_loss')
    train_accuracy = Mean(name='train_accuracy')
    
    # Create a checkpoint object and manager to manage multiple checkpoints
    ckpt = train.Checkpoint(model=training_model, optimizer=optimizer)
    ckpt_manager = train.CheckpointManager(ckpt, "./checkpoints", max_to_keep=3)
    
    # Speeding up the training process
    @function
    def train_step(encoder_input, decoder_input, decoder_output):
        with GradientTape() as tape:
    
            # Run the forward pass of the model to generate a prediction
            prediction = training_model(encoder_input, decoder_input, training=True)
    
            # Compute the training loss
            loss = loss_fcn(decoder_output, prediction)
    
            # Compute the training accuracy
            accuracy = accuracy_fcn(decoder_output, prediction)
    
        # Retrieve gradients of the trainable variables with respect to the training loss
        gradients = tape.gradient(loss, training_model.trainable_weights)
    
        # Update the values of the trainable variables by gradient descent
        optimizer.apply_gradients(zip(gradients, training_model.trainable_weights))
    
        train_loss(loss)
        train_accuracy(accuracy)
    
    
    for epoch in range(epochs):
    
        train_loss.reset_states()
        train_accuracy.reset_states()
    
        print("\nStart of epoch %d" % (epoch + 1))
    
        start_time = time()
    
        # Iterate over the dataset batches
        for step, (train_batchX, train_batchY) in enumerate(train_dataset):
    
            # Define the encoder and decoder inputs, and the decoder output
            encoder_input = train_batchX[:, 1:]
            decoder_input = train_batchY[:, :-1]
            decoder_output = train_batchY[:, 1:]
    
            train_step(encoder_input, decoder_input, decoder_output)
    
            if step % 50 == 0:
                print(f'Epoch {epoch + 1} Step {step} Loss {train_loss.result():.4f} Accuracy {train_accuracy.result():.4f}')
                # print("Samples so far: %s" % ((step + 1) * batch_size))
    
        # Print epoch number and loss value at the end of every epoch
        print("Epoch %d: Training Loss %.4f, Training Accuracy %.4f" % (epoch + 1, train_loss.result(), train_accuracy.result()))
    
        # Save a checkpoint after every five epochs
        if (epoch + 1) % 5 == 0:
            save_path = ckpt_manager.save()
            print("Saved checkpoint at epoch %d" % (epoch + 1))
    
    print("Total time taken: %.2fs" % (time() - start_time))


if __name__ == '__main__':
    main()


