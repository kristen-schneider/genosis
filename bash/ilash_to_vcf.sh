#!bin/sh

PATH=$PATH:~/plink_linux_x86_64_20220402/

vcf_file="/home/sdp/precision-medicine/data/vcf/ALL.chr8.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf"
python_dir="/home/sdp/precision-medicine/python/scripts/ibd/"
ilash_file="/home/sdp/precision-medicine/data/ibd/iLASH/ALL.chr8.ilash.match"
data_dir="/home/sdp/precision-medicine/data/ibd/iLASH/"
cd $data_dir

i=0
echo "iLASH"
while read -r F1 S1 F2 S2 chr bp_start bp_end snp_s snp_e cm percent;
do

	# make vcf from start and end of cm
	echo "Creating ibd segment vcf..."
	out_vcf=$data_dir"ibd.$i.vcf"
	python $python_dir"ilash_to_vcf.py" --vcf $vcf_file --start $bp_start --end $bp_end --out_vcf $out_vcf

	# run plink on vcf
	echo "Running plink on ibd segment vcf..."
	plink --vcf $out_vcf --genome --out "ibd.$i"
	
	# delete vcf
	echo "Removing ibd segment vcf..."


	i=$((i + 1))

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

