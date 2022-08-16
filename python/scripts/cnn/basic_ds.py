import tensorflow as tf
import sampling

def build_dataset_from_file(CNN_input_file):
    """
    Builds a tensorflow dataset object from file with 3 columns:
    sample1, sample2, IBD_distance
    :param CNN_input_file:
    :return: tensorflow dataset object for input to cnn model
    """
    f = open(CNN_input_file, 'r')
    header = f.readline()
    sample1_data = []
    sample2_data = []
    distances_data = []
    c = 0
    for line in f:
        # print(c)
        currS1 = []
        currS2 = []
        A = line.strip().split()
        for i in A[0]:
            currS1.append(int(i))
        sample1_data.append(currS1)
        for j in A[1]:
            currS2.append(int(j))
        sample2_data.append(currS1)

        distances_data.append(A[2])
        c += 1

    # each set of input data is its own tensor
    sample1_ds = tf.data.Dataset.from_tensor_slices(sample1_data)
    sample2_ds = tf.data.Dataset.from_tensor_slices(sample2_data)
    distances_ds = tf.data.Dataset.from_tensor_slices(distances_data)

    # map distance values to float
    sample1_ds_float = sample1_ds.map(lambda x: int(x))
    sample2_ds_float = sample2_ds.map(lambda x: int(x))
    distances_ds_float = distances_ds.map(lambda x: float(x))

    # zip 3 tensor items together into one dataset
    full_ds = tf.data.Dataset.zip((sample1_ds_float, sample2_ds_float, distances_ds_float))
    # puts 1 element into a batch
    full_ds_batch = full_ds.batch(2)
    return full_ds_batch

def sample_without_replacement(sample_IDs_file, sample_encodings_file,
                               ID_distances_file,
                               encoding_distance_file,
                               num_samples, num_variants):

    print('...getting lists')
    sample_IDs_list = sampling.get_sample_ID_list(sample_IDs_file)
    sample_encodings_list = sampling.get_encoding_list(sample_encodings_file)
    # distances_list = sampling.get_distances_list(pairwise_distances_file)

    ID_index_dict = {ID: s for s, ID in enumerate(sample_IDs_list)}
    ID_encoding_dict = sampling.get_ID_encoding_dict(sample_IDs_list,
                                                     sample_encodings_file)
    print('...making pairwise dict')
    pairwise_dict = sampling.get_pairwise_distances_dict(ID_index_dict,
                                                         ID_distances_file)
    print('...sampling without replacement, converting to tensorflow dataset')
    dataset = tf.data.Dataset.from_generator(
        # the generator has to be callable and take no args
        # that's annoying, but we can wrap the call to sample_pairs in a lambda
        lambda: sampling.sample_pairs(sample_encodings_list, pairwise_dict),
        # dtypes of the tensors -- need to adapt to real data
        (tf.int32, tf.int32, tf.float32),
        # shapes of tensors -- need to adapt to real data
        (tf.TensorShape((num_variants,)), tf.TensorShape((num_variants,)), tf.TensorShape(())),
    ).batch(2)

    # for s1, s2, d in dataset:
    #     print(f"{s1.numpy() = }, {s2.numpy() = }, {d.numpy() = }")

    return dataset
