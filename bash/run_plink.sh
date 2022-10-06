#!bin/sh

PATH=$PATH:~/plink_linux_x86_64_20220402/

data_dir="/home/sdp/precision-medicine/data/ped_sim_data/segments/"
cd $data_dir

for vcf in $data_dir*.vcf;
do
	echo $vcf
	filename=$(basename $vcf)
	seg_name=${filename%.*}
	echo $seg_name
	#sed '10q;d' $vcf
	plink --vcf $vcf --genome --out $seg_name
done

