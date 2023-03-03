#!bin/sh

#PATH=$PATH:~/plink_linux_x86_64_20220402/

data_dir="/home/sdp/pmed-local/data/SAS/SAS_sim/"
vcf="/home/sdp/pmed-local/data/SAS/SAS_sim/SAS_pois.vcf.gz"
out="SAS_pois"

plink2 --vcf $vcf \
	--make-king-table \
	--out $data_dir"plink-kingship/"$out
