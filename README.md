# 1. MAMBA ENVIRONMENT
### creating environment from environment.yml (recommended)
`mamba env create -f environment.yml`
### creating raw environment with mamba
`mamba create -c conda-forge -c bioconda -n precision-medicine snakemake htslib plink faiss tensorflow`<br>
<br>
# 2. CONFIG FILE
- edit cpp/configs/segment\_config<br>
- edit config.yaml<br>
<br>
# 3. SNAKEMAKE
run `snakemake`

