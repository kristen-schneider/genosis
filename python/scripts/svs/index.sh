echo "1. building svs index slices..."
start_index=$(date +%s.%3N)

python build_svs_index.py \
	--emb_dir /home/sdp/pmed-local/data/chrm10-12_embeddings/ \
	--idx_dir /home/sdp/pmed-local/data/chrm10-12_svs_index/ \
	--db_samples /home/sdp/pmed-local/data/1kG_hap_IDs.txt \
	--emb_ext emb

end_index=$(date +%s.%3N)
index_time=$(echo "scale=3; $end_index - $start_index" | bc)
echo "--SVS INDEX: $index_time seconds"
