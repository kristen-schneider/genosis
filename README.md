## ABOUT [TO DECIDE NAME]
[ TO DECIDE NAME ] is a tool to query over biobank-scale genotype data and report KNN for a query individual (or set of query individuals).<br>
### Installation
TO UPDATE (for now, see below)<br>
### [Helper Scripts](https://github.com/kristen-schneider/precision-medicine/blob/main/helper_scripts/)
In our paper [ ], and at some points in this README, we refer to helper scripts which were part of our research process and could be useful to you while using this tool.<br>

### Test on [example data](https://github.com/kristen-schneider/precision-medicine/blob/main/example/)
##### 0. Clone the latest version of the repository and update submodules OR pull the latest changes.
```
git clone https://github.com/kristen-schneider/precision-medicine.git
cd precision-medicine
git submodule init
git submodule update
```
```
cd precision-medicine
git pull
``` 
##### 1. Create and activate the mamba environment.
```
mamba env create -f environment.yml
mamba activate pmed
```
##### 2. Run the provided example script with toy data.
I have added a very small set of example data in the directory `data/example_data/`.<br>
This should be an easy way to test that software is working before running on your own dataset..<br>
The snakemake file should already be fit to use example data.<br>
```
snakemake -c1
```

### Run on real data

#### 0. Prepare input data.
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
- A single VCF file for all samples and chromosomes (can be less than full genome, if smaller regions are necessary). See [`example/example.vcf`](https://github.com/kristen-schneider/precision-medicine/blob/main/example/example.vcf.gz) for reference.<br.
**Prepare input MAP file.** `./my_dir/my_map.txt`:
. `./my_dir/data.map`: map file (format: `chromosome` `variant identifier` `position in centimorgans` `base-pair coordinate`). This file will serve as reference to the interpolate map script; which will create a new map (`interpolated.map`) that is used for down stream work.
- see [`helper_scripts/reorder_map.sh`](https://github.com/kristen-schneider/precision-medicine/blob/main/helper_scripts/reorder_map.sh) to fix formatting errors with your map file.
- see [`example/example.map`](https://github.com/kristen-schneider/precision-medicine/blob/main/example/examplemap), or visit [this resource](https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/) for ready-to-download map files for GRCh36, 37, and 38.
4. **TODO** (later development for list of query and database individuals)

#### 1. Edit yaml file to include correct file and directory paths.
- see [`config_snakemake.yaml`](https://github.com/kristen-schneider/precision-medicine/blob/main/example/config_snakemake.yaml) for example.
- We suggest you only modify the following options to reflect appropriate file paths (all other options should be okay to keep as they are).
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
#### 2. Edit snakemake file to point to new configs.
```
configfile: "./my_dir/config_snakemake.yaml"
```

#### 3. Running with snakemake
```
snakemake -c1
```
