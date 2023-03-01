data_dir="/home/sdp/pmed-local/data/"
ped_sim_dir="/home/sdp/pmed-local/ped-sim/"

$ped_sim_dir"ped-sim" \
	-d $data_dir"SAS/SAS_sim/SAS.def" \
	-i $data_dir"SAS/SAS.vcf.gz" \
	-m $data_dir"chr8.map" \
   	-o $data_dir"SAS/SAS_sim/SAS_pois" \
	--pois \
	--founder_ids \
	--keep_phase \
	--fam \
	--bp \
	--mrca
	#--retain_extra -1 \
	#--intf $data_dir"chr8.interfere.tsv" \
