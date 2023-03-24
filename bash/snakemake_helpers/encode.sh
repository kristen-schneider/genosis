#!/bin/bash

bin=$1
config=$2
seg_dir=$3

for vcf_seg in $seg_dir*.vcf; do
	filename=$(basename $vcf_seg)
	seg_name=${filename%.*}
	echo $seg_name
	$bin $config $vcf_seg $seg_dir$seg_name.gt $seg_dir$seg_name.pos $seg_dir$seg_name.af &
done
touch $seg_dir"encode.log"
