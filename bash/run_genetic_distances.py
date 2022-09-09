#!/bin/sh

python_dir='/home/sdp/precision-medicine/python/scripts/pysam'
map_file='/home/sdp/precision-medicine/data/maps/ALL.chr8.interpolated.map'

snp_segment_sizes=(1000 5000 10000 50000 100000 1000000)
cm_segment_sizes=(1 2 3)

out_dir='/home/sdp/precision-medicine/data/snps'
chr=8


echo "Measuring centimorgan lengths for varying segment sizes (SNPs)."
for i in "${snp_segment_sizes[@]}"
do
	echo "...segment size = $i SNPs"
	out_file="$out_dir/chr$chr.snp_size.$i.cM"
	echo "...writing to: $out_file"
	python $python_dir/"count_seg_bp.py" --map $map_file --snps $i > $out_file
done

echo "Counting number of SNPs in each segment for vary segment sizes (cM)."
for i in "${cm_segment_sizes[@]}"
do
        echo "...segment size = $i cMs"
        out_file="$out_dir/chr$chr.cm_size.$i.cM"
        echo "...writing to: $out_file"
        python $python_dir/"count_cm_snps.py" --map $map_file --cm_max $i > $out_file
done
