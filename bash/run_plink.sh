#!bin/sh

source ~/miniconda3/etc/profile.d/conda.sh
conda activate precision-medicine
PATH=$PATH:/plink_linux_x86_64_20220402/

data_dir="/home/sdp/precision-medicine/data/segments/seg_5000/"
cd $data_dir

for vcf in $data_dir*.vcf;
do
	echo $vcf
	plink --vcf $vcf --genome
	filename=$(basename $vcf)
	seg_name=${filename%.*}
	mv $data_dir"plink.genome" $data_dir$seg_name".genome"
done

