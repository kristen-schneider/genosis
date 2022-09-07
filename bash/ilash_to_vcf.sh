#!bin/sh

PATH=$PATH:~/plink_linux_x86_64_20220402/

python_dir="/home/sdp/precision-medicine/python/scripts/ibd/"
ilash_file="/home/sdp/precision-medicine/data/ibd/iLASH/ALL.chr8.match"
data_dir="/home/sdp/precision-medicine/data/ibd/iLASH/"
cd $data_dir

i=0
while IFS="" read -r p || [ -n "$p" ]
do
	start=$(awk '{print $6}')
	end=$(awk '{print $7}')
	#i=$((i + 1))
	echo $end
	#python write_vcf.py --start $start --end $end --out_vcf

done < $ilash_file


#for vcf in $data_dir*.vcf;
#do
#	echo $vcf
#	filename=$(basename $vcf)
#	seg_name=${filename%.*}
#	echo $seg_name
#	sed '10q;d' $vcf
#	#plink --vcf $vcf --genome --out $seg_name
#done

