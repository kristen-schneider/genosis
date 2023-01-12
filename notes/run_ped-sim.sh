data_dir="/home/sdp/pmed-local/data/"
ped_sim_dir="/home/sdp/pmed-local/ped-sim/"

$ped_sim_dir"/ped-sim" \
	-d def_files/EUR.def \
	-i $data_dir"1kGP_high_coverage_Illumina.chr8.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
	-m $data_dir"plink/plink.chr8.GRCh38.map" \
   	-o $data_dir"EUR/EUR_output" \
	--intf $ped_sim_dir"interfere/nu_p_campbell.tsv"
