#!bin/sh

PATH=$PATH:~/plink_linux_x86_64_20220402/

python_dir="/home/sdp/precision-medicine/python/scripts/ibd/"
ilash_file="/home/sdp/precision-medicine/data/ibd/iLASH/ALL.chr8.ilash.match"
test_file="/home/sdp/precision-medicine/data/ibd/iLASH/test.test"
data_dir="/home/sdp/precision-medicine/data/ibd/iLASH/"
cd $data_dir

i=0
while read -r F1 S1 F2 S2 chr bp_start bp_end snp_s snp_e cm percent;
do
	echo "start: " $bp_start
	echo $bp_end
done < $ilash_file

#while read line; do
#	echo $line
	#start=$(awk '{print $6 $7}')
	#echo $start
	#i=$((i + 1))
	#python write_vcf.py --start $start --end $end --out_vcf

#done < $ilash_file


#for vcf in $data_dir*.vcf;
#do
#	echo $vcf
#	filename=$(basename $vcf)
#	seg_name=${filename%.*}
#	echo $seg_name
#	sed '10q;d' $vcf
#	#plink --vcf $vcf --genome --out $seg_name
#done

