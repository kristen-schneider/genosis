import argparse
import sys

import tensorflow as tf
from datasets import PairsDataset


def evaluate_model(args: argparse.Namespace):
    model = tf.keras.models.load_model(args.model_path)
    samples = [line.rstrip().split("\t") for line in open(args.test_samples, "r")]
    test_data = PairsDataset(
        keep_samples=samples,
        sample_id_filename=args.samples_ids,
        sample_pair_filename=args.sample_pairs,
        genotype_filename=args.genotypes,
        shuffle=False,
        batch_size=args.batch_size,
    )
    model.compile(run_eagerly=True)
    model.evaluate(test_data, verbose=1)


if __name__ == "__main__":
    model_path = sys.argv[1]
    test_samples = sys.argv[2]  # 1 sample per line

    parser = argparse.ArgumentParser()
    parser.add_argument("--model_path", type=str, dest="model_path", default=model_path)
    parser.add_argument("--test_samples", type=str, dest="test_samples")
    parser.add_argument("--sample_ids", type=str, dest="samples_ids")
    parser.add_argument("--sample_pairs", type=str, dest="sample_pairs")
    parser.add_argument("--genotypes", type=str, dest="genotypes")
    parser.add_argument("--batch_size", type=int, default=32, dest="batch_size")
    args = parser.parse_args()
    evaluate_model(args)
