## Test on [example data](https://github.com/kristen-schneider/precision-medicine/blob/main/example/)
### Download the latest release and update submodules.
[download latest release](https://github.com/kristen-schneider/precision-medicine/releases/tag/v1.3)
```
cd precision-medicine
git submodule init
git submodule update
```
### Create and activate the mamba environment.
```
mamba env create -f environment.yml
mamba activate pmed
```
###  Run the provided example script with toy data.
I have added a very small set of example data in the directory [`data/example_data/`](https://github.com/kristen-schneider/precision-medicine/blob/main/example/config_snakemake.yaml).<br>
This should be an easy way to test that software is working before running on your own dataset.<br>
Make changs to conda directory (and any other necessary paths) in the config file.<br>
Mostly, the snakemake file should already be fit to use example data.<br>
```
snakemake -c1
```
____________________________________________

## Run on real data
### Prepare input data.
**Create a working directory** `my_dir` with the following structure:
```
my_dir/
	↳log/
	↳segments/
	↳my_vcf.vcf.gz
	↳my_vcf.vcf.gz.tbi
	↳my_map.txt
	↳config_snakemake.yaml
```
**Prepare input VCF file.** `./my_dir/my_vcf.vcf.gz`:
- A single VCF file for all samples and chromosomes. VCF file should be bgzipped and tabixed. Can be less than full genome (e.g. 1 chromosome), if smaller regions are being searched. See [`example/example_merge.vcf.gz`](https://github.com/kristen-schneider/precision-medicine/blob/main/example/example_merge.vcf.gz) for an example VCF.<br>

**Prepare input MAP file.** `./my_dir/my_map.map`:
- map file format (column ordering): `chromosome` `variant identifier` `position in centimorgans` `base-pair coordinate`.
- This file will serve as reference to the interpolate map script; which will create a new map (`interpolated.map`) that is used for down stream work.
- see [`helper_scripts/reorder_map.sh`](https://github.com/kristen-schneider/precision-medicine/blob/main/helper_scripts/reorder_map.sh) to fix formatting errors with your map file.
- see [`example/example_merge.map`](https://github.com/kristen-schneider/precision-medicine/blob/main/example/examplemap), or visit [this resource](https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/) for ready-to-download map files for GRCh36, 37, and 38.

### Edit yaml file to include correct file and directory paths.
- see [`config_snakemake.yaml`](https://github.com/kristen-schneider/precision-medicine/blob/main/example/config_snakemake.yaml) for example.
- We suggest you only modify the following options to reflect appropriate file paths (all other options should be okay to keep as they are, but double checking never hurts 🧐).
 ```
conda_dir:
  /home/sdp/miniconda3/envs/pmed/

data_dir:
  example/
vcf_file:
  example/vcfs/example.8-10.vcf.gz
ref_map:
  example/maps/example.8-10.map
out_dir:
  example/segments/

k:
  20
delim:
  space
```
### Edit snakemake file to point to new configs.
```
configfile: "./my_dir/config_snakemake.yaml"
```

### Running with snakemake
```
snakemake -c1
```
<br>
<br>

____________________________________________
TODO: ADD BEFORE PUBLICATION
- **About section.** Name of tool. Describtion of tool *is a tool to query over biobank-scale genotype data and report KNN for a query individual (or set of query individuals)*.<br>
- **Installation section.** <br>
- **[Helper Scripts](https://github.com/kristen-schneider/precision-medicine/blob/main/helper_scripts/)**.
... helper scripts which were part of our research process and could be useful to you while using this tool.<br>
- **Indexing step. Query step**
- **Images and figures.**
- **Footnotes/acknowledgments.**
- **[link to publication]().**
- **citation for publication.**
____________________________________________

