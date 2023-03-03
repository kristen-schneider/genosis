#!/bin/sh

# path to directories
src_dir="/home/sdp/precision-medicine/cpp/src/"
include_cpp="/home/sdp/precision-medicine/cpp/include/"
include_conda="/home/sdp/miniconda3/envs/pmed/include/"
lib_conda="/home/sdp/miniconda3/envs/pmed/lib/"
include_htslib="/home/sdp/precision-medicine/lib/htslib/"
bin_dir="/home/sdp/precision-medicine/cpp/bin/"
config_file="/home/sdp/precision-medicine/example/config_ex_snakemake.yaml"

bin=$bin_dir"slice-one"

# compile
echo "...compiling"
g++ $src_dir"main_slice.cpp" \
	$src_dir"slice_vcf.cpp" \
	$src_dir"read_map.cpp" \
	$src_dir"utils.cpp" \
	$src_dir"read_config.cpp" \
	-I $include_cpp \
	-I $include_htslb \
	-lhts \
	-o $bin

# execute
echo "...executing."
$bin $config_file
