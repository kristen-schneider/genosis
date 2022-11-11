#!/bin/sh

python_dir="/home/sdp/precision-medicine/python/scripts/plotting"
data_dir="/home/sdp/precision-medicine/data/segments/chr8/seg_5000/"
base_name="chr8.seg."
samples_dir="/home/sdp/precision-medicine/data/samples"
plink_ext=".genome"
faiss_enc_ext=".enc_faissL2"
faiss_emb_ext=".emb_faissL2"
out_dir="/home/sdp/precision-medicine/python/scripts/plotting/png/"
num_seg=50

python $python_dir/distributions.py \
	--data_dir $data_dir \
	--base_name $base_name \
	--plink_ext $plink_ext \
	--faiss_enc_ext $faiss_enc_ext \
	--faiss_emb_ext $faiss_emb_ext \
	--out_dir $out_dir \
	--num_seg $num_seg \
	--train $samples_dir"/train_samples.txt" \
	--test $samples_dir"/test_samples.txt"
