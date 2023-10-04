#!/usr/bin/env bash

input_IDs=$1

while IFS= read -r sample_ID;
do
	echo $sample_ID"_0"
	echo $sample_ID"_1"
done < $input_IDs
