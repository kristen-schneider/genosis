import argparse
import plot_segment_information as psi
import matplotlib.pytplot as plt

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--boundary_file', type=str, required=True, help='Segment boundary file')
    parser.add_argument('--encoding_dir', type=str, required=True, help='Directory containing encodings')
    parser.add_argument('--embedding_dir', type=str, required=True, help='Directory containing embeddings')
    parser.add_argument('--png', type=str, required=True, help='Path to save the plot')

def main():
    args = parse_args()
    boundary_file = args.boundary_file
    encoding_dir = args.encoding_dir
    embedding_dir = args.embedding_dir
    png = args.png

    # 1. count the number of segments for each chromosome
    chrm_segment_count = psi.count_cm(boundary_file)
    # 2. count the length of genotype encodings for each segment
    gt_segment_lengths = psi.count_vector_lengths(encoding_dir, '.gt')
    # 3. count the length of positional encodings for each segment
    pos_segment_lengths = psi.count_vector_lengths(encoding_dir, '.pos')
    # 4. count the length of embeddings for each segment
    embedding_segment_lengths = psi.count_vector_lengths(embedding_dir, '.emb')

    # plot the results
    psi.plot_dimensionality_reduction(gt_segment_lengths,
                                      pos_segment_lengths,
                                      embedding_segment_lengths,
                                      png)