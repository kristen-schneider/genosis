#!/usr/bin/env bash

set -e pipefail

pmed_dir="/home/name/precision-medicine/"
vcf_file="/home/name/precision-medicine/example/example_merge.vcf.gz"
out_dir=$pmed_dir"example/"

# go to project directory
cd $pmed_dir

# load conda and activate snakemake env for run
. /opt/conda/etc/profile.d/conda.sh
conda activate plink

# plink genome
plink \
	--vcf $vcf_file \
	--genome \
	--out $out_dir"plink-genome"

# plink2 kingship
plink2 \
	--vcf $vcf_file \
	--make-king-table
