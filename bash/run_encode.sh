#!/bin/sh

# path to directories
src_dir="/home/sdp/precision-medicine/src/"
include_dir="/home/sdp/precision-medicine/include/"
conda_dir="/home/sdp/miniconda3/envs/faiss/"
bin_dir="/home/sdp/precision-medicine/bin/"
data_dir="/home/sdp/precision-medicine/data/"

# path to vcf and encoded file
vcf_file=$data_dir"vcf/ten.vcf"
encoded_file=$data_dir"encoded/new.encoded.txt"

# search and encoding info
numVariants=9   # number of variants in encoded file
numSamples=2548   # number of samples in encoded file

#source ~/miniconda3/etc/profile.d/conda.sh 
#conda activate faiss
#export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:$conda_dir"lib"

echo "Encoding file: " $vcf_file

bin=$bin_dir"encode"

g++ $src_dir"main_encode.cpp" \
	$src_dir"readVCF.cpp" \
	$src_dir"utils.cpp" \
	-I $include_dir \
	-I $conda_dir"include/" \
	-L $conda_dir"lib" \
	-lhts \
	-o $bin

$bin $vcf_file $encoded_file
