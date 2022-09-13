#!/bin/sh

# path to directories
src_dir="/home/sdp/precision-medicine/cpp/src/"
include_dir="/home/sdp/precision-medicine/cpp/include/"
conda_dir="/home/sdp/miniconda3/envs/pm/"
bin_dir="/home/sdp/precision-medicine/cpp/bin/"
data_dir="/home/sdp/precision-medicine/data/"

# path to vcf and encoded file
config_file="/home/sdp/precision-medicine/cpp/configs/sample.config"

echo "Slicing file: " $vcf_file

bin=$bin_dir"test_slice"

g++ $src_dir"slice_vcf.cpp" \
	$src_dir"read_config.cpp" \
	-I $include_dir \
	-I $conda_dir"include/" \
	-L $conda_dir"lib" \
	-lhts \
	-o $bin

$bin $config_file
