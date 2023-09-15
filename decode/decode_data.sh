#!/usr/bin/env bash

set -e pipefail

### TODO: 
### modify these paths

out_dir="/path/to/out_dir/"

samples="/path/to/sample_list.pns"
chi="/path/to/chi.chi"
out_vcf=$out_dir"out.vcf"

cm_py="/path/to/addCM2map.py"
full_map="/path/to/full_map.map"
cmap="/path/to/cmap.cmap.gz"
out_map=$out_dir"out.map"

###
###

# make VCF file
chitools view -p $samples $chi -f vcf > $out_vcf
bgzip $out_vcf
tabix -p vcf $out_vcf".gz"

# make map file
awk '{ print $1 $3 $4 $2}'
python $cm_py $full_map $cmap
