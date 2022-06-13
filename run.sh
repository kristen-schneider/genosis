#!/bin/sh

# path to directories
src_dir="/home/sdp/precision-medicine/src/"
include_dir="/home/sdp/precision-medicine/include/"
conda_dir="/home/sdp/miniconda3/envs/faiss/"
bin_dir="/home/sdp/precision-medicine/bin/"
data_dir="/home/sdp/precision-medicine/data/"

# path to encoded and query file
encoded_file=$data_dir"encoded/short.encoded.txt"
quireies_file=$data_dir"queries/short.queries.txt"

# search and encoding info
numVariants=9   # number of variants in encoded file
numSamples=3;   # number of samples in encoded file
numQueries=1;   # number of queries in queries file
k=9;            # number of nearest neighbors to report
segmentLengthStart=3;
segmentLengthEnd=9;
segmentLengthStep=1;

source ~/miniconda3/etc/profile.d/conda.sh 
conda activate faiss
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:$conda_dir"lib"
echo "Starting All Experiments."
# for 3 differnt kinds of indexes
for index in {0..0}
do
    	echo "INDEX: " $index

    	#for segment lengths 100-1000 (100 step)
	for segment_length in {3..9..3}
	#for segment_length in {$segmentLengthStart..$segmentLengthEnd}
    	do
        	echo "SEGMENT LENGTH: " $segment_length

        	bin=$bin_dir"index"$index"-length"$segment_length

        	g++ $src_dir"main_faiss.cpp" \
			$src_dir"buildIndex.cpp" \
			$src_dir"searchIndex.cpp" \
			$src_dir"readEncoding.cpp" \
			-I $include_dir \
			-I $conda_dir"include/" \
			-L $conda_dir"lib" \
			-lfaiss \
			-o $bin
        	$bin

    	done
   	echo
done
