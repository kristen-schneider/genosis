#!/usr/bin/env bash

#SBATCH -p short
#SBATCH --job-name=svs-index
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=64gb
#SBATCH --time=1:00:00
#SBATCH --mail-type=NONE
#SBATCH --mail-user=krsc0813@colorado.edu
#SBATCH --output=/Users/krsc0813/precision-medicine/python/scripts/svs/svs-index.out
#SBATCH --error=/Users/krsc0813/precision-medicine/python/scripts/svs/svs-index.err

echo "1. building svs index slices..."
start_index=$(date +%s.%3N)

emb_dir="/Users/krsc0813/precision-medicine/example/embeddings/"
idx_dir="/Users/krsc0813/precision-medicine/example/svs_index/"
ids="/Users/krsc0813/precision-medicine/example/query_hap_IDs.txt"

python build_svs_index.py \
	--emb_dir $emb_dir \
	--idx_dir $idx_dir \
	--db_samples $ids \
	--emb_ext emb

end_index=$(date +%s.%3N)
index_time=$(echo "scale=3; $end_index - $start_index" | bc)
echo "--SVS INDEX: $index_time seconds"
