#!/bin/sh

data_dir='/home/sdp/pmed-local/data/1KG/segments/'

echo "seg SNPs"
for vcf_f in $data_dir*.vcf
do
	filename=$(basename -- $vcf_f)
        base=${filename%.*}
	seg_i=${base##*.}
	#echo $vcf_f
	start_bp=$(sed -n '112p' $vcf_f | awk '{print $2}')
	end_bp=$(tail -n 1 $vcf_f | awk '{print $2}')
	echo $seg_i $start_bp $end_bp
	#SNPs=$(wc -l < $vcf_f)
	#echo $seg_i " " $SNPs
done
