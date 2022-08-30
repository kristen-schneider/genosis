#!/bin/sh

python_dir="/home/sdp/precision-medicine/python/scripts/plotting"
data_dir="/home/sdp/precision-medicine/data/segments/chr0_test/seg_5000/"
plink="genome"
faiss_enc="enc_faissL2"
faiss_emb="emb_faissL2"
k=20
train_samples="/home/sdp/precision-medicine/data/samples/train_samples.txt"
test_samples="/home/sdp/precision-medicine/data/samples/test_samples.txt"
out_dir="/home/sdp/precision-medicine/python/scripts/plotting/"
num_seg=50

python $python_dir/write_stats.py --data_dir $data_dir --plink $plink --faiss_enc $faiss_enc --faiss_emb $faiss_emb --k $k --train $train_samples --test $test_samples --out_dir $out_dir --num_seg $num_seg
