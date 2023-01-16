### RUN ON TEST DATA

##### 0. Clone the latest version of the repository and update submodules.
```
git clone https://github.com/kristen-schneider/precision-medicine.git
cd precision-medicine
git submodule init
git submodule update
```

##### 1. Create and activate the mamba environment.
```
mamba env create -f environment.yml
mamba activate pmed
```

##### 2. Run the provided example script with toy data.
I have added a very small set of example data in the directory `data/example_data/`.<br>
This should be an easy way to test that software is working without having to make changes yet.<br>
The snakemake file should already be fit to use example data.<br>
test with:<br>
```
snakemake -c1
```

### RUN ON REAL DATA

##### 0. Prepare input data.
###### Generate encodings (part 1)
- create a new directory to house your data `mkdir ./my_dir`
- `./my_dir/data.vcf`: single chromosome vcf file for all samples. see [`example/example.vcf`](https://github.com/kristen-schneider/precision-medicine/blob/main/example/example.vcf).
- `./my_dir/data.map`: map file (format: chr start cm end). see [`example/example.map`](https://github.com/kristen-schneider/precision-medicine/blob/main/example/examplemap), or visit [links at bottom of page](https://github.com/kristen-schneider/precision-medicine/edit/main/README.md#more-information-about-data-included-in-config-files) for more map file resources.
- `./my_dir/samples_query.txt`: list of sample IDs against which to query. see [`example/samples_query.txt`](https://github.com/kristen-schneider/precision-medicine/blob/main/example/samples_query.txt) for a single query; add rows for multiple query.
- create a new directory to house your segments `mkdir ./my_dir/segments`

##### 1. Edit yaml files for snakemake options.
###### [`config_ex_snakemake.yaml`](https://github.com/kristen-schneider/precision-medicine/blob/main/example/config_ex_snakemake.yaml)
- copy `config_ex_snakemake.yaml` into your new directory: `cp example/config_ex_snakemake.yaml ./my_dir/`
- Modify the following options to reflect appropriate file paths.
  ```
  home_dir:
    /path/to/project/directory/
  conda_dir:
    /path/to/miniconda3/envs/pm/

  data_dir:
    /path/to/my_dir/
  vcf_file:
    /path/to/my_dir/file.vcf
  query_file:
    /path/to/my_dir/samples_query.txt
  out_dir:
    /path/to/my_dir/segments/
  config_file
    /path/to/my_dir/config_new.yaml
  ```
###### [`config_ex.yaml`](https://github.com/kristen-schneider/precision-medicine/blob/main/example/config_ex.yaml)
- copy `config_ex.yaml` into your new directory: `cp example/config_ex.yaml ./my_dir/`
- Modify the following options to reflect appropriate file paths.
  ```
  home_dir:
  ./
  conda_dir:
  ~/miniconda3/envs/pm
  .
  .
  .
  data_dir:
  example/
  vcf_file:
  example/example.vcf
  map_file:
  example/example.map
  encoding_file:
  example/example_encoding.txt
  sample_IDs_file:
  example/samples_IDs.txt
  query_file:
  example/samples_query.txt
  out_dir:
  example/segments/
  out_base_name:
  example.data
  segment_size:
  1
  ```
##### 2. Edit snakemake file to point to new configs.
```
configfile: "./my_dir/config_new_snakemake.yaml"
```

##### 3. Running with snakemake
```
snakemake -c1
```

### More information about data included in config files.
`map_file`: a genetic map for the data (contains centimorgan distance information). [HapMap genetic maps in cM units from Beagle](https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh38.map.zip)<br> Can use [ilash_analyzer](https://github.com/roohy/ilash_analyzer/blob/master/interpolate_maps.py) to help create map files.<br>
`out_dir`: a subdirectory of `data_dir` which will contain all intermediate segment files.<br>
`out_base_name`: used for naming scheme of all intermediate segment files.

