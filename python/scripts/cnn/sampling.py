# Needs Python >= 3.9
import random
from collections import defaultdict
from typing import Any, Iterable, Mapping, Sequence

import tensorflow as tf

# copied partially from itertools recipes
def grouper(iterable: Iterable, n: int, in_order=True) -> Iterable[tuple[Any, ...]]:
    """
    Collect data into non-overlapping fixed-length chunks or blocks.
    eg: grouper('ABCDEFG', 3) --> ABC DEF

    Note: if n doesn't divide evenly into the length of the iterable,
    the remainder is ommited.

    If in_order is True, return the group in sorted order (ascending).

    """
    if in_order:
        iterable = sorted(iterable)
    args = [iter(iterable)] * n
    return zip(*args)


def sample_without_replacement(n: int) -> list[int]:
    """
    Sample n random numbers without replacement.
    """
    return random.sample(range(n), n)


def sample_pairs(samples: Sequence, distances: Mapping[int, Mapping[int, float]]):
    """
    Yield random pairs of samples and distances
    - samples: like a list or np array containing the samples
    - distances: two level mapping of sample indices to pairwise distance.
    distances[i][j] is the distance between sample i and sample j.

    Note: I am assuming that i < j.

    Also we are assuming that samples fits into memory.
    """
    for i, j in grouper(sample_without_replacement(len(samples)), 2):
        yield samples[i], samples[j], distances[i][j]


def get_ID_encoding_dict(sample_IDs, sample_encodings_file):
    ID_encoding_dict = dict.fromkeys(sample_IDs)

    ID_i = 0
    f = open(sample_encodings_file, 'r')
    for line in f:
        encoding_i = line.strip()
        ID_encoding_dict[sample_IDs[ID_i]] = encoding_i
        ID_i += 1
    f.close()

    if len(sample_IDs) != len(ID_encoding_dict.keys()):
        print('ERROR: not the same number of sample IDs and encodings.')
    return ID_encoding_dict

def get_sample_ID_list(sample_ID_file):
    all_sample_IDs = []
    f = open(sample_ID_file, 'r')
    for line in f:
        ID = line.strip()
        all_sample_IDs.append(ID)
    f.close()

    return all_sample_IDs

def get_encoding_list(sample_encodings_file):
    all_sample_encodings = []
    f = open(sample_encodings_file, 'r')
    for line in f:
        encoding = line.strip()
        encoding_list_ints = [int(i) for i in encoding]
        all_sample_encodings.append(encoding_list_ints)
    f.close()

    return all_sample_encodings

def get_distances_list(pairwise_distances_file):
    distances_list = []
    header = None
    f = open(pairwise_distances_file, 'r')
    for line in f:
        if header == None: header = line
        else:
            A = line.strip().split()
            distance = float(A[2])
            distances_list.append(distance)
    f.close()
    return distances_list

def get_pairwise_distances_dict(ID_index_dict, ID_distances_file):
    pairwise_dict = defaultdict(dict)
    # for k in pairwise_dict.keys():
    #     pairwise_dict[k] = dict()

    f = open(ID_distances_file, 'r')
    header = None
    for line in f:
        if header == None:
            header = line
        else:
            A = line.strip().split()
            sample1_ID = A[0]
            sample1_i = ID_index_dict[sample1_ID]

            sample2_ID = A[1]
            sample2_i = ID_index_dict[sample2_ID]

            distance = float(A[2])
            pairwise_dict[sample1_i][sample2_i] = distance
    f.close()

    return pairwise_dict


# # ------------------------------------------------------------------------------
# # Test it out
# # ------------------------------------------------------------------------------
# if __name__ == "__main__":
#     # i < j ensures that we dont get redundant pairs i, j and j, i.
#     dummy_encodings = ['0p0', '1p0', '2p0', '3p0', '4p0',
#                        '5p0', '6p0', '7p0', '8p0', '9p0']
#     distances = defaultdict(dict)
#     for i in range(10):
#         for j in range(10):
#             if i < j:
#                 distances[i][j] = random.random()
#
#     # Sanity check: 10 choose 2 = 45, so there should be 45 possible pairs
#     # Sampling from 10 pairs of items without replacement should give us 5 pairs.
#     # n = 0
#     # for i in distances:
#     #     for j in distances[i]:
#     #         print(f"{n}: {i=}, {j=}, {distances[i][j]=}")
#     #         n += 1
#
#     dataset = tf.data.Dataset.from_generator(
#         # the generator has to be callable and take no args
#         # that's annoying, but we can wrap the call to sample_pairs in a lambda
#         lambda: sample_pairs(dummy_encodings, distances),
#         # dtypes of the tensors -- need to adapt to real data
#         (tf.string, tf.string, tf.float32),
#         # shapes of tensors -- need to adapt to real data
#         (tf.TensorShape([]), tf.TensorShape([]), tf.TensorShape([])),
#     )
#
#     for s1, s2, d in dataset:
#         print(f"{s1.numpy() = }, {s2.numpy() = }, {d.numpy() = }")
