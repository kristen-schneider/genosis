#!/bin/sh

python_dir='/home/sdp/precision-medicine/python/scripts/distance/'
data_dir='/home/sdp/precision-medicine/data/ped_sim_data/GBR/segments/'

sampleID='fam1_g3-b1-i1'
haps=(0 1)

for encoded_f in $data_dir*.encoded
do
	filename=$(basename -- $encoded_f)
        base=${filename%.*}
	echo "Computing distances for " $filename
	
	for h in "${haps[@]}"
	do
		echo "...haplotype" $h
        	kdist_f=$data_dir$base".hap."$h".kdist"
		edist_f=$data_dir$base".hap."$h".edist"
		svdist_f=$data_dir$base".hap."$h".svdist"
		rdp_dist_f=$data_dir$base".hap."$h".rdpdist"
		python $python_dir"compute_trios.py" \
			--encoded_file $encoded_f \
			--hap $h \
			--query $sampleID \
			> $edist_f
	done
done

