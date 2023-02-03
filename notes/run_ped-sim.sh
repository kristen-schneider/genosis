data_dir="/home/sdp/pmed-local/data/"
ped_sim_dir="/home/sdp/pmed-local/ped-sim/"

$ped_sim_dir"ped-sim" \
	-d $data_dir"EUR/EUR_small.def" \
	-i $data_dir"EUR/EUR_small.vcf.gz" \
	-m $data_dir"maps/chr8.map" \
   	-o $data_dir"EUR/EUR_small_sim" \
	--intf $data_dir"chr8.interfere.tsv" \
	--founder_ids \
	--keep_phase \
	--fam \
	--bp \
	--mrca
