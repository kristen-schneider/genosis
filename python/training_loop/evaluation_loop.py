import argparse
import sys

import matplotlib.pyplot as plt
import tensorflow as tf
from datasets import PairsDataset
from models import DistanceLayer


def evaluate_model(args: argparse.Namespace):
    model = tf.keras.models.load_model(args.model_path)
    samples = [line.rstrip() for line in open(args.test_samples, "r")]
    print(samples)
    test_data = PairsDataset(
        keep_samples=samples,
        sample_id_filename=args.samples_ids,
        sample_pair_filename=args.sample_pairs,
        genotype_filename=args.genotypes,
        shuffle=False,
        batch_size=args.batch_size,
        repeat=False,
    )
    print(test_data.num_pairs)

    distance = DistanceLayer()
    dpred = []
    dtrue = []
    for i, (s1, s2, d) in enumerate(test_data.ds):
        print(f"{i}/{test_data.num_pairs}")
        for truth, pred in zip(
            d.numpy(),
            tf.keras.layers.Dot(axes=-1, normalize=True)([model(s1), model(s2)]).numpy(),
        ):
            dpred.append(pred)
            dtrue.append(truth)

    plt.scatter(dtrue, dpred)
    plt.xlabel("IBD distance")
    plt.ylabel("Predicted distance")
    plt.show()


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
