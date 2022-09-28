import argparse
import tensorflow as tf

from datasets import SingleDataset

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--sample-list",
        type=str,
        required=True,
        help="Path to the file with the samples to use.",
        dest="sample_list",
    )
    parser.add_argument(
        "--model-path",
        type=str,
        required=True,
        help="Path to the model to use.",
        dest="model_path",
    )
    parser.add_argument(
        "--output-path",
        type=str,
        required=True,
        help="Path to the output file.",
        dest="output_path",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=32,
        help="Batch size to use.",
        dest="batch_size",
    )
    parser.add_argument(
        "--sample_id_filename",
        type=str,
        required=True,
        help="File containing all sample ids. Used to create id to int index mapping.",
        dest="sample_id_filename",
    )
    parser.add_argument(
        "--genotype_filename",
        type=str,
        required=True,
        help="File containing all genotype encoding strings.",
        dest="genotype_filename",
    )
    return parser.parse_args()

def main():
    args = get_args() 
    compute_vectors(args)

def compute_vectors(args):
    """
    Load the model and compute the vectors for the dataset.
    """
    # Load the model
    model = tf.keras.models.load_model(args.model_path)

    # which samples to keep
    keep_samples = [s.rstrip() for s in open(args.sample_list, "r").readlines()]

    # Load the dataset
    dataset = SingleDataset(
        keep_samples=keep_samples,
        sample_id_filename=args.sample_id_filename,
        genotype_filename=args.genotype_filename,
        batch_size=args.batch_size,
        shuffle=False,
        repeat=False,
    )

    # Compute the vectors
    with open(args.output_path, "w") as f:
        for samples, encodings in dataset.ds:
            vectors = model(encodings)
            for s, v in zip(samples, vectors):
                print(s.numpy().decode("utf-8"), *v.numpy(), file=f)

if __name__ == '__main__':
    main()
