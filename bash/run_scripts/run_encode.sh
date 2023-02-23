#!/bin/sh

# path to directories
base_dir="/home/sdp/precision-medicine/"
src_dir="/home/sdp/precision-medicine/cpp/src/"
include_dir="/home/sdp/precision-medicine/cpp/include/"
bin_dir="/home/sdp/precision-medicine/cpp/bin/"
data_dir="/home/sdp/pmed-local/data/1KG/segments/"
htslib_dir="/home/sdp/precision-medicine/lib/htslib/"

# path to vcf and encoded file
config="/home/sdp/pmed-local/data/1KG/config.yaml"
sample_IDs=$data_dir"/"
vcf_file=$data_dir"segment.100.vcf"

bin=$bin_dir"encode"

g++ $src_dir"main_encode.cpp" \
	$src_dir"encode_segment.cpp" \
	$src_dir"read_config.cpp" \
	$src_dir"read_map.cpp" \
	$src_dir"map_encodings.cpp" \
	$src_dir"utils.cpp" \
	-I $include_dir \
	-I $htslib_dir \
	-lhts \
	-o $bin

$bin $config $vcf_file $base_dir"out.gt" $base_dir"out.pos"

#for vcf_f in $data_dir*.vcf; do
#	filename=$(basename -- $vcf_f)
#	base=${filename%.*}
#	encoded_f=$data_dir$base".encoded"
#
#	echo "Encoding file: " $vcf_f
#	$bin $config $sample_IDs $vcf_f $encoded_f
#done

