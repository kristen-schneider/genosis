#! bin/bash
bin_dir="/home/sdp/precision-medicine/cpp/bin/"
vcf="/home/sdp/pmed-local/data/1KG/1kGP_high_coverage_Illumina.chr8.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
map="/home/sdp/pmed-local/data/1KG/chr8.interpolated.map"
queries_file="/home/sdp/pmed-local/data/SAS/SAS_sample_IDs.txt"
out_dir="/home/sdp/precision-medicine/compare/"
iLASH_ex="/home/sdp/iLASH/"
k=20

# ilash
echo "RUNNING iLASH"
ilash_start=`date +%s.%N`
$iLASH_ex"/build/ilash" $iLASH_ex/"1kg-chr8.config" > $out_dir/"1kg-chr8.ilash.out"
ilash_end=`date +%s.%N`
ilash_runtime=$( echo "$ilash_end - $ilash_start" | bc -l )
echo "ilash runtime: $ilash_runtime"
echo 

# hap-ibd to go here

## plink --genome
#echo "RUNNING PLINK --GENOME"
#pg_start=`date +%s.%N`
#plink --vcf $vcf --genome --out $out_dir"plink.1kg-chr8.genome"
#pg_end=`date +%s.%N`
#pg_runtime=$( echo "$pg_end - $pg_start" | bc -l )
#$bin_dir"plink-genome" $out_dir"plink.1kg-chr8.genome.genome" $k $out_dir"topk_genome.SAS" $queries_file
#echo "plink --genome runtime: $pg_runtime"
#
## plink --make-king-table
#echo "RUNNING PLINK --MAKE-KING-TABLE"
#pk_start=`date +%s.%N`
#plink2 --vcf $vcf --make-king-table --out $out_dir"plink.1kg-chr8.king"
#$bin_dir"plink-king" $out_dir"plink.1kg-chr8.king.kin0" $k $out_dir"topk_king.SAS" $queries_file
#pk_end=`date +%s.%N`
#pk_runtime=$( echo "$pk_end - $pk_start" | bc -l )
#echo "plink --make-king runtime: $pk_runtime"

