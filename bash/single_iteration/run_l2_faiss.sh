#!/bin/sh

# to run a single faiss search on a input encoding vectors
# path to directories
src_dir="/home/sdp/precision-medicine/cpp/src/"
include_dir="/home/sdp/precision-medicine/cpp/include/"
conda_dir="/home/sdp/miniconda3/envs/faiss/"
bin_dir="/home/sdp/precision-medicine/cpp/bin/"
data_dir="/home/sdp/pmed-local/data/"

# path to encoded and query file
all_samples=$data_dir"/1KG/all_hap_IDs.txt"
query_samples=$data_dir"/EAS/EAS_hap_IDs.txt"
train_samples=$data_dir"database_IDs.txt"
k=20 # number of nearest neighbors to report
delim="space"

echo "Compiling..."

bin=$bin_dir"faiss_l2"

g++ $src_dir"faiss_l2.cpp" \
	$src_dir"build_index.cpp" \
	$src_dir"read_encodings.cpp" \
	$src_dir"search_index.cpp" \
	$src_dir"faiss_utils.cpp" \
	-I $include_dir \
	-I $conda_dir"include/" \
	-L $conda_dir"lib/" \
	-lfaiss \
	-o $bin

echo "Executing..."
for encoded_f in $data_dir"/1KG/segments/"*.encoded
do
        filename=$(basename -- $encoded_f)
        base=${filename%.*}
        echo "Running FAISS on" $filename
        faiss_out=$data_dir"segments/"$base".norm.faissl2"
	$bin $all_samples $encoded_f $query_samples $encoded_f $k $delim
	#echo $bin $all_samples $encoded_f $all_samples $encoded_f $k $delim > $faiss_out
done


#$bin $train_samples $encoded_file $train_samples $encoded_file $k $delim
#$bin $query_samples $embedded_file $query_samples $embedded_file $k $delim
