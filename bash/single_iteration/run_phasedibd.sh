#!/bin/sh

# path to directories
src_dir="/home/sdp/precision-medicine/python/scripts/ibd/"
conda_dir="/home/sdp/miniconda3/envs/pm/"
data_dir="/home/sdp/precision-medicine/data/"

# path to vcf and hap file
vcf_file=$data_dir"vcf/ALL.chr8.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf"
map_file=$data_dir"maps/ALL.chr8.interpolated.map"
#vcf_file="/home/sdp/phasedibd/tests/data/pedigree_vcf/10.vcf"
#map_file=$data_dir

echo $vcf_file

python $src_dir"phased_ibd.py" --vcf $vcf_file --map $map_file > $data_dir"ibd/phasedibd/ALL.chr8.phasedibd.out"
