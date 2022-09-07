#!bin/sh

PATH=$PATH:~/plink_linux_x86_64_20220402/

vcf_file="/home/sdp/precision-medicine/data/vcf/ALL.chr8.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf"
python_dir="/home/sdp/precision-medicine/python/scripts/ibd/"
phasedibd_file="/home/sdp/precision-medicine/data/ibd/phasedibd/pibd.csv"
data_dir="/home/sdp/precision-medicine/data/ibd/phasedibd/"
cd $data_dir

i=0
echo "phased-ibd"
while read -r chr id1 id2 id1_h id2_h start end start_cm end_cm start_bp end_bp;
do

	# make vcf from start and end of cm
	echo "Creating ibd segment vcf..."
	out_vcf=$data_dir"ibd.$i.vcf"
	python $python_dir"ilash_to_vcf.py" --vcf $vcf_file --start $start_bp --end $end_bp --out_vcf $out_vcf

	# run plink on vcf
	echo "Running plink on ibd segment vcf..."
	plink --vcf $out_vcf --genome --out "ibd.$i"
	
	# delete vcf
	echo "Removing ibd segment vcf..."
	rm $out_vcf
	rm "ibd.$i.log"
	rm "ibd.$i.nosex"


	i=$((i + 1))

done < $phasedibd_file

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

