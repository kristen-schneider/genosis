#!/bin/sh

python_dir="/home/sdp/precision-medicine/python/scripts/plotting"
data_dir="/home/sdp/precision-medicine/data/segments/chr8/seg_5000/"
plink_ext="log"
faiss_enc_ext="enc_faissL2"
faiss_emb_ext="emb_faissL2"
out_dir="/home/sdp/precision-medicine/python/scripts/plotting/txt/"
num_seg=50

python $python_dir/write_times.py \
	--data_dir $data_dir \
	--plink_ext $plink_ext \
	--faiss_enc_ext $faiss_enc_ext \
	--faiss_emb_ext $faiss_emb_ext \
	--out_dir $out_dir \
	--num_seg $num_seg \
