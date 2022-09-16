#!/bin/sh

python_dir='/home/sdp/precision-medicine/python/scripts/distance/'
data_dir='/home/sdp/precision-medicine/data/segments/chr8-30x/'

sampleID='HG00405'
paternalID='HG00403'
maternalID='HG00404'
haps=(0 1)

for encoded_f in $data_dir*.encoded
do
	filename=$(basename -- $encoded_f)
        base=${filename%.*}
	echo "Computing distances for " $filename
	
	for h in "${haps[@]}"
	do
		echo "...haplotype" $h
        	kdist_f=$data_dir$base"hap."$h".kdist"

		python $python_dir"compute_trios.py" \
			--encoded_file $encoded_f \
			--hap $h \
			--query $sampleID \
			> $kdist_f
	done
done

