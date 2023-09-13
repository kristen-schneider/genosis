## SETUP
____________________________________________
**Prepare input VCF file** `data_dir/my_vcf.vcf.gz`:
- A single VCF file for all samples and chromosomes. VCF file should be bgzipped and tabixed. Can be less than full genome (e.g. 1 chromosome), if smaller regions are being searched. See [`example/example_merge.vcf.gz`](https://github.com/kristen-schneider/precision-medicine/blob/main/example/example_merge.vcf.gz) for an example VCF.<br>

**Prepare input MAP file** `data_dir/my_map.map`:
- map file format (column ordering):<br>
`chromosome` `variant identifier` `position in centimorgans` `base-pair coordinate`.
- This file will serve as reference to the interpolate map script; which will create a new map (`interpolated.map`) that is used for down stream work.
- see [`helper_scripts/reorder_map.sh`](https://github.com/kristen-schneider/precision-medicine/blob/main/helper_scripts/reorder_map.sh) to fix formatting errors with your map file.
- see [`example/example_merge.map`](https://github.com/kristen-schneider/precision-medicine/blob/main/example/examplemap), or visit [this resource](https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/) for ready-to-download map files for GRCh36, 37, and 38.
### Download the latest release and update submodules.
```
git clone git@github.com:kristen-schneider/precision-medicine.git
cd precision-medicine
git submodule init
git submodule update
```
### Run with singularity
Build singularity container from recipe file:
```
sudo singularity build pmed.sif pmed_recipe.def
```
Edit file paths:<br>
- Open [example_config.yml](https://github.com/kristen-schneider/precision-medicine/blob/main/example/config_singularity.yml) and change paths appropriately.<br>
- Open [example_run.sh](https://github.com/kristen-schneider/precision-medicine/blob/main/example/config_run.sh) and change path to config file appropriately.<br>
## WORKFLOW
____________________________________________
### Run with singularity _(mount necessary data directories if necessary)_:<br>
```
singularity run --bind /path/to/data_dir/ pmed.sif bash example_run.sh
# if conda is not working properly, try the following #
## $ conda init
## $ . /opt/conda/etc/profile.d/conda.sh
```
