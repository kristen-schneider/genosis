import os 

def count_cm(boundary_file):
    """
    count the number of segments for each chromosome
    """
    # format of boundary file:
    # chr segment start end
    chrm_segment_count = {}
    with open(boundary_file, 'r') as f:
        for line in f:
            line = line.strip().split()
            chrm = line[0]
            if chrm not in chrm_segment_count:
                chrm_segment_count[chrm] = 0
            chrm_segment_count[chrm] += 1
    return chrm_segment_count

def count_vector_lengths(vector_dir, ext):
    """
    count the length of vectors for each segment

    @param vector_dir: directory containing genotype encodings
    @return: dictionary mapping segment id to encoding length
    """
    segment_lengths = {}
    for filename in os.listdir(vector_dir):
        if filename.endswith(ext):
            # filename: chrm1.segment1.ext
            segment_id = filename.replace(ext, '')

            # open the file and count the length of each line
            vector_file = os.path.join(vector_dir, filename)
            all_sample_lengths = count_vector_length(vector_file)
            segment_lengths[segment_id] = all_sample_lengths
    return segment_lengths

def count_vector_length(segment_file):
    """
    count the length of all sample vectors in a segment file

    @param segment_file: file containing genotype encodings for a segment
    @return: list of vector lengths for all samples
    """
    vector_lengths = []
    with open(segment_file, 'r') as f:
        for line in f:
            line = line.strip().split()
            vector_length = len(line) - 1 # first column is the segment id
            vector_lengths.append(vector_length)
    return vector_lengths

def plot_dimensionality_reduction(gt_vector_lengths,
                                  pos_vector_lengths,
                                  emb_vector_lengths,
                                  png):
    """
    plot the dimensionality reduction of genotype encodings to positional encodings to embeddings

    @param gt_vector_lengths: dictionary mapping segment id to genotype encoding lengths
    @param pos_vector_lengths: dictionary mapping segment id to positional encoding lengths
    @param emb_vector_lengths: dictionary mapping segment id to embedding lengths
    @param png: path to save the plot
    """
    
    # x axis is the genotype encoding lengths
    # y axis is the positional and embedding lengths

    # plot the results
    fig, ax = plt.subplots()
    for segment_id in gt_vector_lengths:
        gt_length = gt_vector_lengths[segment_id]
        pos_length = pos_vector_lengths[segment_id]
        emb_length = emb_vector_lengths[segment_id]
        ax.scatter(gt_length, pos_length, marker='o', color='steelblue')
        ax.scatter(gt_length, emb_length, marker='x', color='salmon')
    ax.set_xlabel('Genotype Encoding Length')
    ax.set_ylabel('Vector Length')
    ax.set_title('Dimensionality Reduction')
    ax.legend(['Positional Encoding', 'Embedding'])
    plt.show()
    plt.savefig(png)
