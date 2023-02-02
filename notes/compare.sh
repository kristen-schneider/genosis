#! bin/bash
vcf=""
map=""
out_dir=""

# plink --genome
pg_start=`date +%s.%N`
#plink --vcf $vcf --genome --out $out_dir"plink.genome"
sleep 5s 
pg_end=`date +%s.%N`
pg_runtime=$( echo "$pg_end - $pg_start" | bc -l )
echo "plink --genome runtime: $pg_runtime"

# plink --make-king-table
pk_start=`date +%s.%N`
#plink --vcf $vcf --make-king-table --out $out_dir"plink.king"
sleep 5s 
pk_end=`date +%s.%N`
pk_runtime=$( echo "$pk_end - $pk_start" | bc -l )
echo "plink --make-king runtime: $pk_runtime"

# ilash
ilash_start=`date +%s.%N`
#./build/ilash 1kg-chr8.config > $out_dir"1kg-chr8.out"
sleep 5s 
ilash_end=`date +%s.%N`
ilash_runtime=$( echo "$ilash_end - $ilash_start" | bc -l )
echo "ilash runtime: $ilash_runtime"

# hap-ibd
hapibd_start=`date +%s.%N`
#java -jar hap-ibd.jar gt=$vcf map=$map out=$out
# CURRENT BUG WITH DUPLICATES
sleep 5s 
hapibd_end=`date +%s.%N`
hapibd_runtime=$( echo "$hapibd_end - $hapibd_start" | bc -l )
echo "hapibd runtime: $hapibd_runtime"

