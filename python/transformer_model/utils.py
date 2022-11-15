import numpy as np
import tensorflow as tf
from tensorflow import convert_to_tensor
from typing import Mapping, Tuple

"""
Used to read and split incoming data.
KEY ASSUMPTION: order of genotype encodings is same
    as order of sample IDs.
    TODO: assert this later
"""

def read_file_to_list(filename: str) -> list[str]:
    """
    Reads lines from file and returns all lines in a list.

    :param filename: (str) file name with all data
    :return:
        file_list: (list) list of data from file (str)
    """
    f = open(filename, 'r')
    file_list = f.read().splitlines()
    return file_list

def map_sample_names_to_index(sample_names_file: str) -> Mapping[str, int]:
    """
    A dictionary mapping with key = sample_name, value = index

    :param sample_names_file: (str) file with all sample names
    :return:
        sample_names_index: (map) key = sample name, value = index
    """
    # get list of sample names
    sample_names_list = read_file_to_list(sample_names_file)
    # make dictionary where key is sample name index and value is genotype encoding
    sample_names_index = {sample_name: i*2 for i, sample_name in enumerate(sample_names_list)}
    return sample_names_index

def map_genotype_encoding_to_index(genotype_encoding_file: str) -> Mapping[int, list]:
    """
    A dictionary mapping with key = integer ID, value = genotype encoding

    :param genotype_encoding_file: (str) file with all sample genotype encodings
    :return:
        genotype_encodings_index: (map) key = sample name index,
                                        value = genotype encoding
    """
    # get list of genotyope encodings
    genotype_encodings = [l.strip().split()[1:] for l in read_file_to_list(genotype_encoding_file)]
    # make dictionary where key is index and value is genotype encoding
    genotype_encodings_index = {i: ge for i, ge in
                                enumerate([[np.uint8(j) for j in g]
                                           for g in genotype_encodings])}
    return genotype_encodings_index

def map_positional_encoding_to_index(positional_encoding_file: str) -> Mapping[int, list]:
    """
    A dictionary mapping with key = integer ID, value = positional encoding

    :param positional_encoding_file: (str) file with all sample positional encodings
    :return:
        genotype_encodings_index: (map) key = sample name index,
                                        value = positinoal encoding
    """
    # get list of positional encodings
    positional_encodings = [l.strip().split()[1:] for l in read_file_to_list(positional_encoding_file)]
    # make dictionary where key is index and value is positional encoding
    positional_encodings_index = {i: pe for i, pe in
                                enumerate([[np.uint8(j) for j in p]
                                           for p in positional_encodings])}
    return positional_encodings_index


def split_samples(sample_names_file: Mapping[str, int],
                  train_ratio: float
                  ) -> Tuple[list[str], list[str], list[str]]:
    """
    splits a list of data files into training, testing,
    and validating samples.

    :param genotype_names_index: (str) file with all sample names
    :param train_ratio: percent of samples to be used for training
    :return: training set, testing set, validation set
    """
    # get list of sample names
    sample_names_list = read_file_to_list(sample_names_file)
    num_samples = len(sample_names_list)

    # make training set
    num_training = int(np.floor(train_ratio * num_samples))
    training_sample_names = sample_names_list[:num_training]

    num_remaining = num_samples - num_training
    # make testing set
    num_testing = num_remaining // 2
    testing_sample_names = sample_names_list[num_training: num_training+num_testing]

    # make validating set
    validating_sample_names = sample_names_list[num_training+num_testing:]

    return training_sample_names, testing_sample_names, validating_sample_names

def get_both_haplotypes(sample_names_index: Mapping[str, int] ,
                        input_IDs: list[str]):
    """
    Takes a list of sample IDs and creates a list of indexes
    coresponding to hap0 and hap1 for each sample ID

    :param sample_names_index: dictionary with
                        key: sample name, value: index.
    :param IDs: list of strings which are sample IDs
    """
    sample_idx_hap_0 = [sample_names_index[ID] for ID in input_IDs]
    sample_idx_hap_1 = [sample_names_index[ID]+1 for ID in input_IDs]
    sample_idx_all = sample_idx_hap_0 + sample_idx_hap_1
    return sample_idx_all


def get_tensors(sample_names_index,
                    encodings_index,
                    training_IDs):
    """
    Converts genotype encodings to tensors object.

    :param sample_names_index: dictionary with
                        key: sample name, value: index.
    :param encodings_index: dictionary with
                        key: sample name index, key: encoding
    :param training_IDs: list of sample names for training
    :return: gt_tensors: tensors object of genotypes
    """
    sample_idx_all = get_both_haplotypes(sample_names_index,
                        training_IDs)

    training_encodings = [encodings_index[ID_idx] for ID_idx in sample_idx_all]
    tensors = convert_to_tensor(training_encodings, dtype=tf.int32)
    return tensors
