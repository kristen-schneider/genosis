# SETUP
## Download the latest release and update submodules.
```
git clone git@github.com:kristen-schneider/precision-medicine.git
cd precision-medicine
git submodule init
git submodule update
```
____________________________________________
# RUN EXAMPLE DATA
```
cd precision-medicine
mkdir example/log
mkdir example/err
sbatch example/example_run.sh
```
____________________________________________
# RUN NEW DATA
## Prepare working directory
```
mkdir local_project
cd local_project
mkdir log
mkdir err
```
## Prepare input data
**Prepare input VCF file** `data_dir/my_vcf.vcf.gz`:
- A single VCF file for all samples and chromosomes. VCF file should be bgzipped and tabixed. Can be less than full genome (e.g. 1 chromosome), if smaller regions are being searched. See [`example/example_merge.vcf.gz`](https://github.com/kristen-schneider/precision-medicine/blob/main/example/example_merge.vcf.gz) for an example VCF. <br>
**Prepare input MAP file** `data_dir/my_map.map`:
- map file format (column ordering):<br>
`chromosome` `variant identifier` `position in centimorgans` `base-pair coordinate`.
- This file will serve as reference to the interpolate map script; which will create a new map (`interpolated.map`) that is used for down stream work.
- see [`helper_scripts/reorder_map.sh`](https://github.com/kristen-schneider/precision-medicine/blob/main/helper_scripts/reorder_map.sh) to fix formatting errors with your map file.
- see [`example/example_merge.map`](https://github.com/kristen-schneider/precision-medicine/blob/main/example/examplemap), or visit [this resource](https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/) for ready-to-download map files for GRCh36, 37, and 38.
**Prepare list of sample IDs for database and query:
- One file with a list of sample IDs to include in the database (for building index). See [`example/database_IDs.txt`](https://github.com/kristen-schneider/precision-medicine/blob/main/example/database_IDs.txt) for an example using 1000 Genomes samples. <br>
- One file with a list of sample IDs to include in the query set (for searching index). See [`example/query_IDs.txt`](https://github.com/kristen-schneider/precision-medicine/blob/main/example/query_IDs.txt) for an example using 1000 Genomes samples. <br>

## Make custom run files and edit file paths.
- Make copy of [example_run.sh](https://github.com/kristen-schneider/precision-medicine/blob/main/run/run.sh) and edit paths appropriately.<br>
- Make copy of [example_config.yml](https://github.com/kristen-schneider/precision-medicine/blob/main/example/example.yml) and edit paths appropriately.<br>
- Make copy of [embed_config.yml](https://github.com/kristen-schneider/precision-medicine/blob/main/run/cluster_config.yml) and edit node configurations appropriately.<br>

## Run
```
cd local_project
sbatch local_run.sh
```
____________________________________________

# RUN WITH SINGULARITY

## Build singularity container

Build singularity container from recipe file:
```
singularity build gess.sif gess.def
```
## Run with singularity _(mount necessary data directories if necessary)_
```
singularity run
  --bind /path/to/data_dir/ \
  ./singularity/pmed.sif \
  bash ./run/singularity.sh \
```
____________________________________________

