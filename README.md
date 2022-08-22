## 1. Activating mamba environment
*creating environment from environment.yml (recommended)*<br>
`mamba env create -f environment.yml`<br>
*creating new environment with mamba*<br>
`mamba create -c conda-forge -c bioconda -n precision-medicine snakemake htslib plink faiss tensorflow`<br>
<br>

## 2. Editing config files
- edit `cpp/configs/segment_config`<br>
- edit `config.yaml`<br>
<br>

## 3. Running with snakemake
run `snakemake`

