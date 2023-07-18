echo "...building faiss index slices..."
start_index=$(date +%s.%3N)

emb_dir="/home/sdp/pmed-local/data/chrm10-12_embeddings/"
for emb in $emb_dir*; do
	echo $emb
	python build_faiss_index.py \
		--emb $emb \
		--idx_dir /home/sdp/pmed-local/data/chrm10-12_faiss_index/ \
		--db_samples /home/sdp/precision-medicine/example/database_IDs.txt 
done

end_index=$(date +%s.%3N)
index_time=$(echo "scale=3; $end_index - $start_index" | bc)
echo "--FAISS INDEX: $index_time seconds"
