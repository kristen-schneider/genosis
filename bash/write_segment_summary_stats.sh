#!/bin/sh

python_dir="/home/sdp/precision-medicine/python/scripts/plotting"
data_dir="/home/sdp/precision-medicine/data/segments/chr8/seg_5000/"
plink_ext=".genome"
faiss_enc_ext=".enc_faissHNSW"
faiss_emb_ext=".emb_faissHNSW"
faiss_idx="hnsw"
base_name="chr8.seg."
train_samples="/home/sdp/precision-medicine/data/samples/train_samples.txt"
test_samples="/home/sdp/precision-medicine/data/samples/test_samples.txt"
k=20
out_dir="/home/sdp/precision-medicine/python/scripts/plotting/txt/"
num_seg=50

python $python_dir/write_stats.py \
	--data_dir $data_dir \
	--plink_ext $plink_ext \
	--faiss_enc_ext $faiss_enc_ext \
	--faiss_emb_ext $faiss_emb_ext \
	--base_name $base_name \
	--faiss_idx $faiss_idx \
	--train $train_samples \
	--test $test_samples \
	--k $k \
	--num_seg $num_seg \
	--out_dir $out_dir
