#!/bin/sh

python_dir="/home/sdp/precision-medicine/python/scripts/plotting"
data_dir="/home/sdp/precision-medicine/data/segments/chr8/seg_5000/"
plink_ext="genome"
faiss_enc_ext="enc_faissL2"
faiss_emb_ext="emb_faissL2"
k=20
train_samples="/home/sdp/precision-medicine/data/samples/train_samples.txt"
test_samples="/home/sdp/precision-medicine/data/samples/test_samples.txt"
out_dir="/home/sdp/precision-medicine/python/scripts/plotting/"
num_seg=50
faiss_index='l2'

#python $python_dir/write_stats.py --data_dir $data_dir --plink_ext $plink_ext --faiss_enc_ext $faiss_enc_ext --faiss_emb_ext $faiss_emb_ext --k $k --train $train_samples

python $python_dir/write_stats.py \
	--data_dir $data_dir \
	--plink_ext $plink_ext \
	--faiss_enc_ext $faiss_enc_ext \
	--faiss_emb_ext $faiss_emb_ext \
	--k $k \
	--train $train_samples \
	--test $test_samples \
	--out_dir $out_dir \
	--num_seg $num_seg \
	--faiss_index $faiss_index
