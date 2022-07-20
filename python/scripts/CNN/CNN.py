import utils

from tensorflow import keras
from tensorflow.keras.layers import Conv1D

CNN_input_file = '/Users/kristen/PycharmProjects/genotype-encoding/data/CNN.input.txt'

def main():
    [vector_length, num_vectors] = CNN_input_dimensions(CNN_input_file)

#    CNN(int(vector_length), int(num_vectors))
#
#def get_dimensions(CNN_input_file):
#    vector_length = -1
#    num_vectors = -1
#
#    f = open(CNN_input_file, 'r')
#    header = f.readline()
#    num_vectors = 0
#    for line in f:
#        if num_vectors == 0:
#            d = line.strip().split()
#            sample1_encoded = d[0]
#            sample2_encoded = d[1]
#            if len(sample1_encoded) != len(sample2_encoded):
#                print('Error. Two encodings of differnt length.')
#                return -1
#            else:
#                vector_length = len(sample1_encoded)
#        num_vectors += 1
#    f.close()
#    return vector_length, num_vectors
#
#def CNN(vector_length, num_vectors):
#    model = keras.models.Sequential()
#    model.add(Conv1D(1, kernel_size=vector_length, input_shape=(num_vectors, 1)))
#    model.summary()
#
#if __name__ == '__main__':
#    main()
