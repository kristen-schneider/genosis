#!/bin/sh

python_dir="/home/sdp/precision-medicine/python/scripts/plotting"
data_dir="/home/sdp/precision-medicine/data/segments/chr8/seg_5000/"
base_name="chr8.seg."
plink_ext=".log"
faiss_enc_ext=".enc_faissL2"
faiss_emb_ext=".emb_faissL2"
faiss_idx="l2"
num_seg=50
out_dir="/home/sdp/precision-medicine/python/scripts/plotting/txt/"

python $python_dir/write_times.py \
	--data_dir $data_dir \
	--base_name $base_name \
	--plink_ext $plink_ext \
	--faiss_enc_ext $faiss_enc_ext \
	--faiss_emb_ext $faiss_emb_ext \
	--faiss_idx $faiss_idx \
	--num_seg $num_seg \
	--out_dir $out_dir 
