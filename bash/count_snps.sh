#!/bin/sh

data_dir='/home/sdp/precision-medicine/data/segments/chr8-30x/'

echo "seg SNPs"
for vcf_f in $data_dir*.vcf
do
	filename=$(basename -- $vcf_f)
        base=${filename%.*}
	seg_i=${base##*.}
	SNPs=$(wc -l < $vcf_f)
	echo $seg_i " " $SNPs
done
