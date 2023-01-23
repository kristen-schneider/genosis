import tensorflow as tf
from typing import Mapping, Tuple

import utils

def nchoosek(n: int, k: int) -> int:
    return int(math.factorial(n) / (math.factorial(k) * math.factorial(n - k)))

class PairsDatasetUnsupervised:
    def __init__(
        self,
        #sample_id_filename: str,
        #genotype_filename: str,
        sample_idx_all: list[int],
        genotype_encodings_index: Mapping[int, list],
        training_IDs: list[str],
        shuffle: bool = False,
        encoding_type: tf.DType = tf.int8,
        batch_size: int = 32,
        repeat: bool = False,
    ):
        """
        TF dataset for sample pairs
        no target distances
        """
        #keep_samples = get_sample_names(sample_id_filename)
        #genotypes = load_genotypes_unsupervised(genotype_filename, keep_samples)
        genotypes = [genotype_encodings_index[idx] for idx in sample_idx_all]
        
        self.num_pairs = nchoosek(len(genotypes), 2) // batch_size
        self.num_variants = len(genotypes[0])
        
        '''
        if shuffle:
            generator = lambda: combinations(
                random.sample(genotypes, len(genotypes)), 2
            )
        else:
            generator = lambda: combinations(genotypes, 2)

        self.ds = tf.data.Dataset.from_generator(
            generator,
            (encoding_type, encoding_type),
            (
                tf.TensorShape([self.num_variants]),
                tf.TensorShape([self.num_variants]),
            ),
        )
        if repeat:
            self.ds = self.ds.repeat()
        self.ds = self.ds.batch(batch_size).prefetch(tf.data.AUTOTUNE)
        '''
