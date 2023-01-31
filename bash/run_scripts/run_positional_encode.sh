#!/bin/sh

# path to directories
src_dir="/home/sdp/precision-medicine/cpp/src/"
include_dir="/home/sdp/precision-medicine/cpp/include/"
include_htslib="/home/sdp/precision-medicine/lib/htslib/"
bin_dir="/home/sdp/precision-medicine/cpp/bin/"
data_dir="/home/sdp/pmed-local/data/1KG/segments/"

config_file="/home/sdp/precision-medicine/notes/ancestry_configs/config_1KG.yaml"
vcf_slice="/home/sdp/pmed-local/data/1KG/segments/1KG.data.seg.0.vcf"
pos_encode="test_slice.pos"
#map_file="/home/sdp/pmed-local/data/1KG/chr8.interpolated.map"


bin=$bin_dir"test_pos_encode"
rm $bin
g++ $src_dir"main_positional_encode.cpp" \
	$src_dir"encode_positions.cpp" \
	$src_dir"read_map.cpp" \
	$src_dir"read_config.cpp" \
	$src_dir"map_encodings.cpp" \
	$src_dir"utils.cpp" \
	-I $include_dir \
	-I $include_htslib \
	-lhts \
	-o $bin

$bin $config_file $vcf_slice $pos_encode


#for vcf_f in $data_dir*.vcf; do
#	filename=$(basename -- $vcf_f)
#	base=${filename%.*}
#	encoded_f=$data_dir$base".encoded"
#
#	echo "Encoding file: " $vcf_f
#	$bin $config $sample_IDs $vcf_f $encoded_f
#done

