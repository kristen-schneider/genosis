#!/bin/sh

segments_dir="/home/sdp/pmed-local/data/1KG/segments/"
python_script="/home/sdp/precision-medicine/notes/python-faiss/faiss-indices.py"
k=20

for encoded_f in $segments_dir*".gt"
do
        filename=$(basename -- $encoded_f)
        base=${filename%.*}
        echo "Running FAISS on" $filename
        faiss_out=$segments_dir$base
        python $python_script $segments_dir $encoded_f $base $k
        #$bin $all_samples $encoded_f $query_samples $encoded_f $k $delim
        #echo $bin $all_samples $encoded_f $all_samples $encoded_f $k $delim > $faiss_out
done

