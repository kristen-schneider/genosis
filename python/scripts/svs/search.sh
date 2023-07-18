echo "searching svs index slices..."
start_search=$(date +%s.%3N)

python search_svs_index.py \
	--idx_dir /home/sdp/pmed-local/data/chrm10-12_svs_index/ \
	--emb_dir /home/sdp/pmed-local/data/chrm10-12_embeddings/ \
	--emb_ext emb \
	--db_samples /home/sdp/pmed-local/data/1kG_hap_IDs.txt \
	--q_samples /home/sdp/pmed-local/data/1kG_hap_IDs.txt \
	--k 20 \
	--out_dir /home/sdp/pmed-local/data/chrm10-12_svs_results/

end_search=$(date +%s.%3N)
search_time=$(echo "scale=3; $end_search - $start_search" | bc)
echo "--SVS SEARCH: $search_time seconds"

