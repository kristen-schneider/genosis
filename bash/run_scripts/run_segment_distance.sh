#!/bin/sh

# to run a single faiss search on a input encoding vectors


# path to directories
conda_dir="/home/sdp/miniconda3/envs/pm/"
data_dir="/home/sdp/precision-medicine/data/ped_sim_data/GBR_large/"
python_dir="/home/sdp/precision-medicine/python/scripts/distance/"

# path to encoded and query file
query_samples=$data_dir"test.txt"
encoded_files=$data_dir"segments/"

source ~/miniconda3/etc/profile.d/mamba.sh 
mamba activate pm

export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:$conda_dir"lib"
echo $LD_LIBRARY_PATH

for encoded_f in $encoded_files*.encoded; do
	filename=$(basename $encoded_f);
	seg_name=${filename%.*};
	#echo $encoded_f;
	echo $seg_name;
	dist_file=$encoded_files$seg_name.eucdist;
	#echo $dist_file
	python $python_dir"compute_segment_distance.py" \
		--encoded_file $encoded_f \
		--query_file $query_samples \
		> $dist_file;
done;
