#!/bin/sh

# path to directories
src_dir="/home/sdp/precision-medicine/cpp/src/"
include_dir="/home/sdp/precision-medicine/cpp/include/"
bin_dir="/home/sdp/precision-medicine/cpp/bin/"

# path to config file
config_file="/home/sdp/precision-medicine/cpp/configs/sample.config"

bin=$bin_dir"slice"

# compile
g++ $src_dir"slice_vcf.cpp" \
	$src_dir"read_config.cpp" \
	-I $include_dir \
	-I $conda_dir"include/" \
	-L $conda_dir"lib" \
	-lhts \
	-o $bin

# execute
$bin $config_file
