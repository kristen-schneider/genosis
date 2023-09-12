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
- Run singularity _(mount necessary data directories if necessary)_:<br>
```
singularity run --build /path/to/data_dir/ pmed.sif
conda init
bash example_run.sh
```
