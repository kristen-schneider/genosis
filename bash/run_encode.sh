#!/bin/sh

# path to directories
src_dir="/home/sdp/precision-medicine/cpp/src/"
include_dir="/home/sdp/precision-medicine/cpp/include/"
conda_dir="/home/sdp/miniconda3/envs/pm/"
bin_dir="/home/sdp/precision-medicine/cpp/bin/"
data_dir="/home/sdp/precision-medicine/data/segments/test_slice/"

# path to vcf and encoded file
config="/home/sdp/precision-medicine/cpp/configs/sample.config"
sample_IDs="/home/sdp/precision-medicine/data/samples/chr8-30x.samples"
vcf_file=$data_dir"chr8.seg.0.vcf"
encoded_file=$data_dir"chr8.seg.0.encoded"

echo "Encoding file: " $vcf_file

bin=$bin_dir"test_encode"

g++ $src_dir"encode_vcf.cpp" \
	$src_dir"map_encodings.cpp" \
	$src_dir"read_config.cpp" \
	$src_dir"utils.cpp" \
	-I $include_dir \
	-I $conda_dir"include/" \
	-L $conda_dir"lib" \
	-lhts \
	-o $bin

$bin $config $sample_IDs $vcf_file $encoded_file
