#!/bin/bash

# loading modules
echo "...loading modules."
module load anaconda

# creating and activating conda environment
echo "...activating conda environment."
conda activate pmed

# running encoding segment 
data_dir=""
cpp_dir=""
htslib_dir=""
echo "...encoding fixed length samples."

# compile
# g++
#	cpp/src/main_encode_fixed.cpp 
#	cpp/src/encode_fixed_segment.cpp
#	cpp/src/map_encodings.cpp
#	cpp/src/utils.cpp -I
#	cpp/include/
#	-I lib/htslib/ 
#	-lhts 
#	-o ./fixed

g++ \
	$cpp_dir"src/main_encode_fixed.cpp" \
	$cpp_dir"src/encode_fixed_segment.cpp" \
	$cpp_dir"src/map_encodings.cpp" \
	$cpp_dir"src/utils.cpp" \
	-I $cpp_dir"include/" \
	-I $htslib_dir \
	-lhts \
	-o $cpp_dir"bin/encode-fixed"

# ./fixed 8 test.vcf example/sample_IDs.txt encoding.txt gt fpos
echo ---ENCODING VCF SEGMENTS---
for vcf_f in $data_dir*.vcf.gz; do
	filename=$(basename $vcf_f);
	seg_name=${filename%.vcf.*};
	echo ... $seg_name;
	chrm_idx=${seg_name#*chrm};
	chrm_idx=${chrm_idx%%.*};
	echo ...$chrm_idx;


	#{input.bin} " \
       #        $chrm_idx " \
       #        $vcf_f " \
       #        {config.data_dir}sample_IDs.txt " \
       #        {config.encoding_file} " \
       #        {config.data_dir}interpolated.map " \
       #        {config.out_dir}${{seg_name}}.gt " \
       #        {config.out_dir}${{seg_name}}.pos " \
       #        {config.out_dir}${{seg_name}}.af " \
       #         >> {output.encode_log};" \
done 
