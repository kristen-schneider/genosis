#!/bin/sh

# path to directories
src_dir="/home/sdp/precision-medicine/cpp/src/"
include_dir="/home/sdp/precision-medicine/cpp/include/"
bin_dir="/home/sdp/precision-medicine/cpp/bin/"
data_dir="/home/sdp/precision-medicine/data/ped_sim_data/segments/"

# path to vcf and encoded file
config="/home/sdp/precision-medicine/cpp/configs/sample.config"
sample_IDs="/home/sdp/precision-medicine/data/ped_sim_data/simulated.samples"
#vcf_file=$data_dir"chr8-30x.seg.9.vcf"
#encoded_file=$data_dir"chr8-30x.seg.9.encoded"

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

for vcf_f in $data_dir*.vcf; do
	filename=$(basename -- $vcf_f)
	base=${filename%.*}
	encoded_f=$data_dir$base".encoded"

	echo "Encoding file: " $vcf_f
	$bin $config $sample_IDs $vcf_f $encoded_f
done

