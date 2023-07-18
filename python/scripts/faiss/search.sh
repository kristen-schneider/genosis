echo "searching faiss index slices..."
start_search=$(date +%s.%3N)

idx_dir="/home/sdp/pmed-local/data/chrm10-12_faiss_index/"
for idx in $idx_dir*; do
        echo $idx
	python search_faiss_index.py \
		--idx $idx \
		--emb_dir /home/sdp/pmed-local/data/chrm10-12_embeddings/ \
		--emb_ext emb \
		--k 20 \
		--database_samples /home/sdp/precision-medicine/example/database_IDs.txt \
		--query_samples /home/sdp/precision-medicine/example/database_IDs.txt \
		--out_dir /home/sdp/pmed-local/data/chrm10-12_faiss_results/

done


end_search=$(date +%s.%3N)
search_time=$(echo "scale=3; $end_search - $start_search" | bc)
echo "--FAISS SEARCH: $search_time seconds"

