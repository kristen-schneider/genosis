#! bin/bash
vcf=""
map=""
out_dir=""

# plink --genome
plink --vcf $vcf --genome --out $out_dir"plink.genome"

# plink --make-king-table
plink --vcf $vcf --make-king-table --out $out_dir"plink.king"

# ilash
./build/ilash 1kg-chr8.config > $out_dir"1kg-chr8.out"

# hap-ibd
java -jar hap-ibd.jar gt=$vcf map=$map out=$out
