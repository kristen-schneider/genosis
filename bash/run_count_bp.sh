#!/bin/sh

python_script='/home/sdp/precision-medicine/python/scripts/pysam/count_seg_bp.py'
map_file='/home/sdp/precision-medicine/data/maps/ALL.chr8.interpolated.map'
segment_sizes=(1000 5000 10000 50000 100000 1000000)
out_dir='/home/sdp/precision-medicine/data/snps'
chr=8

for i in "${segment_sizes[@]}"
do
	echo "counting cM for $i SNPs"
	out_file="$out_dir/chr$chr.seg_size.$i.cM"
	echo "...writing to: $out_file"
	python $python_script --map $map_file --snps $i > $out_file
done
