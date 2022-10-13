TEST_SAMPLES="/home/sdp/precision-medicine/data/samples/chr8-30x/testing.samples"
PYTHON_DIR="/home/sdp/precision-medicine/python/scripts/distance/"
DATA_DIR="/home/sdp/precision-medicine/data/segments/chr8-30x/"
ENCODING_FILE="/home/sdp/precision-medicine/data/segments/chr8-30x/chr8-30x.seg.86.encoded"
HAPS=(0 1)

OUT_FILE=$DATA_DIR"chr8-30x.seg.86.test.edist"
rm $OUT_FILE
while read sample; do
	echo "$sample"
	echo "$sample" >> $OUT_FILE
	for h in "${HAPS[@]}"
	do
		echo $h
		python $PYTHON_DIR"compute_trios.py" \
                        --encoded_file $ENCODING_FILE \
                        --hap $h \
                        --query $sample \
                        >> $OUT_FILE
	done

done < $TEST_SAMPLES
