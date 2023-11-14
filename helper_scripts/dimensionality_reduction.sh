#!/usr/bin/env bash

#SBATCH --partition=short
#SBATCH --job-name=dimensionality-reduction
#SBATCH --output=/Users/krsc0813/1KG_data/chr1_22/out/dim-red.out
#SBATCH --error=/Users/krsc0813/1KG_data/chr1_22/err/dim-red.err
#SBATCH --time=0-23:00:00
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mail-type=NONE
#SBATCH --mail-user=krsc0813@colorado.edu

pmed_dir="/Users/krsc0813/precision-medicine/"
cd $pmed_dir/python/scripts/supplemental/

python plot_segment_information.py \
    --boundary_file "/Users/krsc0813/1KG_data/chr1_22/segment_boundary.map"\
    --encoding_dir "/Users/krsc0813/1KG_data/chr1_22/encodings/" \
    --embedding_dir "/Users/krsc0813/1KG_data/chr1_22/embeddings/" \
    --png "/Users/krsc0813/1KG_data/chr1_22/segments.png"
