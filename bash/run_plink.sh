#!bin/sh

PATH=$PATH:~/plink_linux_x86_64_20220402/

data_dir="/home/sdp/precision-medicine/data/ibd/iLASH/"
cd $data_dir

for vcf in $data_dir*.vcf;
do
	echo $vcf
	filename=$(basename $vcf)
	seg_name=${filename%.*}
	echo $seg_name
	#plink --vcf $vcf --genome --out $seg_name
done

