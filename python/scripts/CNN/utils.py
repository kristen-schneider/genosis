def CNN_input_dimensions(tf_record):
    vector_length = -1
    num_vectors = -1

    f = open(tf_record, 'r')
    header = f.readline()
    num_vectors = 0
    for line in f:
        if num_vectors == 0:
            d = line.strip().split()
            sample1_encoded = d[0]
            sample2_encoded = d[1]
            if len(sample1_encoded) != len(sample2_encoded):
                print('Error. Two encodings of differnt length.')
                return -1
            else:
                vector_length = len(sample1_encoded)
        num_vectors += 1
    f.close()
    return vector_length, num_vectors

