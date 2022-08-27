"""
Used to load data from sample files.
KEY ASSUMPTION: The order of samples is the same in the sample IDs file and
the genotype encodings file.
"""
import random
from pprint import pprint
from typing import Container, Iterable, Mapping, Tuple

import numpy as np
import tensorflow as tf


def get_sample_names(filename: str) -> list[str]:
    """
    Loads the sample names from file
    """
    with open(filename, "r") as f:
        return f.read().splitlines()


def load_sample_ids(filename: str) -> Mapping[str, int]:
    """
    Loads the sample IDs from a file.
    :param filename: The file to load the sample IDs from.
    :return: A dict mapping sample IDs to indices.
    """
    with open(filename, "r") as f:
        sample_ids = f.read().splitlines()
    return {sample_id: i for i, sample_id in enumerate(sample_ids)}


# TODO figure out using memory mapped file to load this
def load_genotypes(filename: str) -> list[list[np.uint8]]:
    """
    load file with genotype encodings and returns a dict keyed by int id
    """
    with open(filename, "r") as f:
        genotypes = f.read().splitlines()  # string of ints with no delimiter
    return [[np.uint8(i) for i in list(g.split()[1])] for g in genotypes]


# get all pairs of samples IDs along with their target distance
def get_sample_pairs(
    filename: str,
) -> list[Tuple[str, str, float]]:
    """
    Loads the sample pairs from a file.
    :param sample_ids: The sample IDs to use.
    :param filename: The file to load the sample pairs from.
    :return: A dict mapping sample pairs to target distances.
    """
    sample_pairs: list[Tuple[str, str, float]] = []

    with open(filename, "r") as f:
        next(f)  # skip header
        for line in f:
            s1, s2, d = line.split()
            sample_pairs.append((s1, s2, float(d)))
    return sample_pairs


# for each sample, get all pairs and yeild their genotype vectors and target distance
def yield_sample_pair(
    sample_ids: Mapping[str, int],
    sample_pairs: list[Tuple[str, str, float]],
    genotypes: list[list[np.uint8]],
    shuffle: bool = False,
) -> Iterable[tuple[list[np.uint8], list[np.uint8], float]]:
    """
    Yields the sample pairs and their target distances.
    :param sample_ids: sample IDs -> integer IDs.
    :param target_dict: sample pair -> distance.
    :param genotypes: integer IDs -> genotype vectors.
    :return: A generator yielding the sample pairs and their target distances.
    """
    if shuffle:
        random.shuffle(sample_pairs)
    for s1, s2, d in sample_pairs:
        yield (genotypes[sample_ids[s1]], genotypes[sample_ids[s2]], d)


def partition_samples(
    train_ratio: float,
    sample_file: str,
) -> Tuple[list[str], list[str], list[str]]:
    """
    Partitions the samples into training, validation and test sets.
    :param train_ratio: The ratio of samples to use for training.
    :param samples: The samples to partition.
    :return: A tuple of lists of sample IDs for each set.
    """
    samples = get_sample_names(sample_file)
    n = len(samples)

    # take the firsn train_ratio% of samples for training
    train_n = np.floor(train_ratio * len(samples)).astype(int)
    train_samples = samples[:train_n]

    # the rest is split 50-50 into val and test
    leftover_n = n - train_n
    val_n = leftover_n // 2
    val_samples = samples[train_n : train_n + val_n]
    test_samples = samples[train_n + val_n :]

    return train_samples, val_samples, test_samples


# TODO generator cant take advantage of multiprocessing due to race conditions
# figure this out later if the GPU doesn't get data fast enough.
# TODO add option for mapping to onehot.
class PairsDataset:
    def __init__(
        self,
        keep_samples: Container,  # leave out pairs that have samples outside this set
        sample_id_filename: str,
        sample_pair_filename: str,
        genotype_filename: str,
        shuffle: bool = False,
        encoding_type: tf.DType = tf.int8,
        target_type: tf.DType = tf.float32,
        batch_size: int = 32,
        repeat: bool = True,
    ):
        """
        TF dataset for the sample pairs and target distances.
        """
        sample_ids = load_sample_ids(sample_id_filename)
        genotypes = load_genotypes(genotype_filename)
        sample_pairs = get_sample_pairs(sample_pair_filename)

        # filter out pairs that have samples outside the keep_samples set
        sample_pairs = [
            (s1, s2, d)
            for s1, s2, d in sample_pairs
            if s1 in keep_samples and s2 in keep_samples
        ]

        self.num_pairs = len(sample_pairs) // batch_size
        self.num_variants = len(genotypes[0])

        self.ds = tf.data.Dataset.from_generator(
            lambda: yield_sample_pair(sample_ids, sample_pairs, genotypes, shuffle),
            (encoding_type, encoding_type, target_type),
            (
                tf.TensorShape([self.num_variants]),
                tf.TensorShape([self.num_variants]),
                tf.TensorShape([]),
            ),
        )
        if repeat:
            self.ds = self.ds.repeat()
        self.ds = self.ds.batch(batch_size).prefetch(10)


class SingleDataset:
    def __init__(
        self,
        keep_samples: Iterable,  # leave out if samples outside this set
        sample_id_filename: str,
        genotype_filename: str,
        shuffle: bool = False,
        # encoding_type: tf.DType = tf.int8,
        batch_size: int = 32,
        repeat: bool = True,
    ):
        """
        TF dataset for single samples.
        Each element is a tuple of (sample_id, genotype_vector)
        """
        # map sample IDs to indices
        sample_ids = load_sample_ids(sample_id_filename)

        # integer indexed list of genotype encodings
        genotypes = load_genotypes(genotype_filename)

        # filter out genotype strings that have samples outside the keep_samples set
        genotypes = [genotypes[sample_ids[s]] for s in keep_samples]

        # create dataset that yeilds tuples of (sample_id, genotype_vector)
        self.ds = tf.data.Dataset.from_tensor_slices((keep_samples, genotypes))

        self.num_variants = len(genotypes[0])
        # self.ds = tf.data.Dataset.from_tensor_slices(data)
        if repeat:
            self.ds = self.ds.repeat()
        if shuffle:
            self.ds = self.ds.shuffle(buffer_size=len(genotypes))
        self.ds = self.ds.batch(batch_size).prefetch(10)


if __name__ == "__main__":
    # sample_ids = load_sample_ids("/home/murad/data/toy_model_data/ALL.sampleIDs")
    # genotypes = load_genotypes("/home/murad/data/toy_model_data/chr0.seg.0.encoding")
    # sample_pairs = get_sample_pairs("/home/murad/data/toy_model_data/chr0.seg.0.cnn")
    # pprint(sample_pairs)
    S = SingleDataset(
        keep_samples=[
            "HG00122",
            "HG00123",
            "HG00125",
        ],
        sample_id_filename="/home/murad/data/toy_model_data/ALL.sampleIDs",
        genotype_filename="/home/murad/data/toy_model_data/chr0.seg.0.encoding",
        shuffle=False,
        batch_size=32,
        repeat=False,
    )
    for x in S.ds:
        print(x.shape)
        exit()

