# DATA
### VCF FILES
-[1KG 30x-grch38](https://www.internationalgenome.org/data-portal/data-collection/30x-grch38)<br>
-[phased VCFs](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/)<br>
-[sample information](https://www.internationalgenome.org/api/beta/sample/_search/igsr-1000%20genomes%2030x%20on%20grch38.tsv.tsv)<br>
-[pedigree information](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/1kGP.3202_samples.pedigree_info.txt)<br>
-vcf to vcf.gz to vcf.gz.tbi
```
bcftools view -S subset_sample_IDs.txt full.vcf.gz > subset_sample.vcf
bgzip -c subset_sample.vcf > subset_sample.vcf.gz
tabix -p vcf subset_sample.vcf.gz
```
-[PED-SIM, simulated](https://github.com/williamslab/ped-sim)
- following instructions on the README.md, all necessary dependencies should be installed with mamba.<br>
- ped-sim Makefile needs a simple change to include mamba environment include directory. copy of this edit is [here](https://github.com/kristen-schneider/precision-medicine/tree/main/notes/Makefile).<br>
- def files for subpopulations are included [here](https://github.com/kristen-schneider/precision-medicine/tree/main/notes/def_files)
- bash file to run a ped-sim simulation [here](https://github.com/kristen-schneider/precision-medicine/tree/main/notes/def_files)
- bash file to run a ped-sim simulation [here](https://github.com/kristen-schneider/precision-medicine/tree/main/notes/run_ped-sim.sh)
### MAP FILES
- [HapMap genetic maps in cM units from Beagle](https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh38.map.zip)<br>
- [ilash analyzer](https://github.com/roohy/ilash_analyzer/blob/master/interpolate_maps.py) to help create map files.<br>
```
# create ped/map files from VCF using plink
plink --vcf /path/to/vcf --recode 01 --output-missing-genotype . --out
# interpolate map file with genetic distance using HapMap or OMNI map
python python/scripts/utils/interploate_map.py --map /path/to/map --ref_map /path/to/ref_map  --out_map /path/to/out_map
```

# KING
```
plink2 --vcf /path/to/vcf.vcf.gz --make-king-table`
```

# CREATE VCF SUBSET FOR TESTING
```
# write header to smaller file
bcftools view -h full_vcf.vcf > small_vcf.vcf 
# write first set of variants to file
tabix full_vcf.vcf chr1:0-2000000 >> small_vcf.vcf
# bgzip and tabix vcf file
bgzip small_vcf.vcf
tabix -p vcf small_vcf.vcf
```

# MERGE VCF FILES
```
bcftools concat vcf1.vcf.gz vcf2.vcf.gz vcf3.vcf.gz > out.vcf
```

# iLASH MAP
```
paste <(awk '{print $1 $2}' 1kG-chrm8.map) <(awk '{print $2}' interpolated.map) <(awk '{print $4}' 1kG-chrm8.map) > ilash-chrm8.map
```

# CONDA/MAMBA 
```
mamba create -n pmed -c bioconda -c conda-forge -c pytorch -c defaults bcftools boost bzip2 faiss-cpu gxx htslib make numpy plink plink2 pysam python snakemake vcftools
```
### CHANNELS
- bioconda
- conda-forge
- defaults
### DEPENDENCIES
- bcftools
- boost
- bzip2
- faiss-cpu
- gxx
- htslib
- make
- numpy
- plink
- pysam
- python
- snakemake
- vcftools
