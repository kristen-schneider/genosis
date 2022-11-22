#!/bin/sh

# normalizes all output encodings 


# path to directories
conda_dir="/home/sdp/miniconda3/envs/pm/"
data_dir="/home/sdp/precision-medicine/data/ped_sim_data/GBR_large/"
python_dir="/home/sdp/precision-medicine/python/scripts/utils/"

# path to encoded and query file
query_samples=$data_dir"self.txt"
encoded_files=$data_dir"segments/"

source ~/miniconda3/etc/profile.d/mamba.sh 
mamba activate pm

export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:$conda_dir"lib"
echo $LD_LIBRARY_PATH

cd $python_dir
#python "normalize_encodings.py" --encoded_file $encoded_files"GBR_large.simulated.seg.0.encoded" > 'test.txt'

for encoded_f in $encoded_files*.encoded; do
	filename=$(basename $encoded_f);
	seg_name=${filename%.*};
	echo $seg_name;
	norm_file=$encoded_files$seg_name.norm;
	#echo $dist_file
	python $python_dir"normalize_encodings.py" \
		--encoded_file $encoded_f \
		> $norm_file;
done;
