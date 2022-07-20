import utils

from tensorflow import keras
from tensorflow.keras.layers import Conv1D

CNN_input_file = '/Users/kristen/PycharmProjects/genotype-encoding/data/CNN.input.txt'

def main():
    [vector_length, num_vectors] = CNN_input_dimensions(CNN_input_file)

#    CNN(int(vector_length), int(num_vectors))
#
#def CNN(vector_length, num_vectors):
#    model = keras.models.Sequential()
#    model.add(Conv1D(1, kernel_size=vector_length, input_shape=(num_vectors, 1)))
#    model.summary()
#
#if __name__ == '__main__':
#    main()
