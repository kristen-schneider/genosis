#!/bin/sh

# path to directories
src_dir="/home/sdp/precision-medicine/src/"
include_dir="/home/sdp/precision-medicine/include/"
conda_dir="/home/sdp/miniconda3/envs/faiss/"
bin_dir="/home/sdp/precision-medicine/bin/"
data_dir="/home/sdp/precision-medicine/data/"

# path to encoded and query file
encoded_file=$data_dir"encoded/test.encoded.txt"
queries_file=$data_dir"queries/test.queries.txt"

# search and encoding info
numVariants=9
#2548903   # number of variants in encoded file
numSamples=15
#2548   # number of samples in encoded file
numQueries=1   # number of queries in queries file
k=15
#2548;            # number of nearest neighbors to report
segmentLength=5
#500

source ~/miniconda3/etc/profile.d/conda.sh 
conda activate faiss
#export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:$conda_dir"lib/"
#echo $LD_LIBRARY_PATH

bin=$bin_dir"one"

g++ $src_dir"main_faiss.cpp" \
	$src_dir"buildIndex.cpp" \
	$src_dir"searchIndex.cpp" \
	$src_dir"readEncoding.cpp" \
	-I $include_dir \
	-I $conda_dir"include/" \
	-L $conda_dir"lib/" \
	-lfaiss \
	-o $bin 

$bin $encoded_file $queries_file $numVariants $numSamples $numQueries $k $segmentLength

