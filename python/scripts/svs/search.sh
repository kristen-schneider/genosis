#!/usr/bin/env bash

#SBATCH -p short
#SBATCH --job-name=svs-search
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=64gb
#SBATCH --time=1:00:00
#SBATCH --mail-type=NONE
#SBATCH --mail-user=krsc0813@colorado.edu
#SBATCH --output=/Users/krsc0813/precision-medicine/python/scripts/svs/svs-search.out
#SBATCH --error=/Users/krsc0813/precision-medicine/python/scripts/svs/svs-search.err


echo "searching svs index slices..."
start_search=$(date +%s.%3N)

emb_dir="/Users/krsc0813/precision-medicine/example/embeddings/"
idx_dir="/Users/krsc0813/precision-medicine/example/svs_index/"
res_dir="/Users/krsc0813/precision-medicine/example/svs_results/"
ids="/Users/krsc0813/precision-medicine/example/query_hap_IDs.txt"

python search_svs_index.py \
	--idx_dir $idx_dir \
	--emb_dir $emb_dir \
	--emb_ext emb \
	--db_samples $ids \
	--q_samples $ids \
	--k 20 \
	--out_dir $res_dir

end_search=$(date +%s.%3N)
search_time=$(echo "scale=3; $end_search - $start_search" | bc)
echo "--SVS SEARCH: $search_time seconds"

