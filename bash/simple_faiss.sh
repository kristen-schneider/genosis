#!/bin/sh

# to run a single faiss search on a input encoding vectors


# path to directories
src_dir="/home/sdp/precision-medicine/cpp/src/"
include_dir="/home/sdp/precision-medicine/cpp/include/"
conda_dir="/home/sdp/miniconda3/envs/precision-medicine/"
bin_dir="/home/sdp/precision-medicine/cpp/bin/"
data_dir="/home/sdp/precision-medicine/data/"

# path to encoded and query file
encoded_file="/home/sdp/genotype-encoding/data/segments/seg_1000/ALL.chr14.seg.0.encoding"
queries_file="/home/sdp/genotype-encoding/data/segments/seg_1000/ALL.chr14.seg.0.encoding"
#encoded_file=$data_dir"encoded/short.encoded.txt "
#queries_file=$data_dir"queries/short.queries.txt "
k=2548 # number of nearest neighbors to report

source ~/miniconda3/etc/profile.d/conda.sh 
conda activate precision-medicine
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:$conda_dir"lib"

echo "Compiling..."
bin=$bin_dir"simple_faiss"
g++ $src_dir"single_faiss.cpp" \
	$src_dir"buildIndex.cpp" \
	$src_dir"readEncoding.cpp" \
	$src_dir"searchIndex.cpp" \
	$src_dir"utils.cpp" \
	-I $include_dir \
	-I $conda_dir"include/" \
	-L $conda_dir"lib" \
	-lfaiss \
	-o $bin

echo "Executing..."

$bin $encoded_file $queries_file $k
