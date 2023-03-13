#!/bin/bash

bin=$1
hap_IDs=$2
seg_dir=$3

for enc_seg in $seg_dir*.gt; do
	filename=$(basename $enc_seg)
	seg_name=${filename%.*}
	echo $seg_name
	$bin $hap_IDs $enc_seg $seg_dir$seg_name.faiss.idx &
done
touch $seg_dir"build_index.log"
