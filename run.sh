#!/bin/sh

src_dir="/home/sdp/precision-medicine/src/"
bin_dir="/home/sdp/precision-medicine/bin/"
data_dir="/home/sdp/precision-medicine/data/"
encoded_file=$data_dir"encoded/short.encoded.txt"


echo "Starting All Experiments."

# for 3 differnt kinds of indexes
for index in {0..3}
do
    echo "INDEX: " $index

    # for segment lengths 100-1000 (100 step)
    for segment_length in {100..1000..100}
    do
        echo "SEGMENT LENGTH: $segment_length"

        bin=$bin_dir"index"$index"-length"$segment_length

        g++ $src_dir"main.cpp" -o $bin  # compile code
        $bin $index $segment_length     # execute binary

    done
    echo
done
