
# Running with snakemake

## Generating training data and training a model
**TODO** need to modify this to work with slurm.

The pipeline started out as data exploration, but ended up becoming
an all in one training set generation/model training pipeline. Located in
`exploration_pipeline.smk`, with configuration at `conf/exp_config.yaml`.

#### Config parameters
##### Data prep config
- `random_seed`: fixes the random number generator for better reproducibility during.
  I've set it to 42.
- `segment_prefix`, `gt_ext`, `pos_ext`, `file_segment_delimiter`, `seg_offset`: These are
  here to allow parseing of filenames that have a more arbitrary format just incase they change
  from the current `segment.$N.{pos|gt}`.
- `num_pairings`: n pairings for each sample from each subpopulation during sampling
- `samples_list`: This is the list of all samples (`${sample_name}_${haplotype}`).
   One way to do it is to take one of the `.gt` files and do `cut -d' ' -f1 $gt_file > list.txt`
- `sample_table`: This is the 1kg super/sub populations table. (not currently in use,
   but maybe someday)
- `segments_dir`: Directory containing the `.gt` files
- `outdir`: output directory of all the results

##### Model training config
- **TODO** add options for half precision training, etc.
- `model_prefix`: used to name runs in `wandb` and also dictates the name of the directory
   containing model checkpoints (`$model_prefix.checkpoints`).
- `train_method`: How to train model.  Currently just implemented siamese.
- `model_type`: Model architecture for backbone.  Currently just conv1d.
- `loss_fn`: loss function for training. `mse` seems best for now.
- `batch_size`: My workstation gpu (12GB VRAM) dictates that we use at most 32 for now,
   higher mem GPUs might allow for more.
- `grad_accum`: To achieve a higher effective batch size without running out of memory,
  we can run multiple batches and accumulate the gradients before making the backprop update.
  This trades memory for time.  For example when `batch_size` is 32 and `grad_accum` is 16,
  the effective batch size is 512.
- `n_workers`: number of processes for the training data loader.  Set this as high as possible
  so that the gpu never has to wait for a batch to be served up.
- `n_epochs`: maximum number of epochs to train for.  I set it at 100.
- `early_stop_patience`: number of validation checks without improvement to wait for before
  stopping training before it has reached the max number of epochs.
- `lr`: base learning rate of the optimizer. Currently set to 1e-3, but can experiment
  with marginally higher or lower values.
- `weight_decay`: weight decay regularization factor.  Higher will be stronger regularization.
  Use if the model seems to overfit quickly.  Doesn't seem to be the case for our current small
  model, but could happen if we increase the number of parameters of the model.
- `n_layers`: number of layers deep the model is.
- `dropout`: dropout probability for dropout layers in the model.  Haven't made use of this
  yet due to the fact we aren't overfitting yet, but use sparingly since it doesn't play
  well with batch norm layers.
- `stride`: *Convnets only* - stride of convolution operation.  Haven't messed with this yet.
- `kernel_size`: *Convnets only* - size of convolutional kernels.  A way to increase the model's receptive field.
- `padding`: *Convnets only* - wether or not to pad vectors with zeros or not to make the
  output the same size as the input.  For padding use `same`.  Otherwise use `valid`.

## Model evaluation
The eval pipeline is located in `eval_pipeline.smk` with coniguration at
`conf/eval_config.yaml`.

#### Config parameters
- `model`: path to a trained model checkpoint (`.ckpt`)
- `gpu`: Set this to true if you want to run inference on GPU.  It takes only a few
  seconds to encode the test set with gpu, but the pipeline also supports cpu parallelization
  so it won't be super slow without it.  Just use more cores in the snakemake call.
- `samples_list`: This is the list of all samples (`${sample_name}_${haplotype}`).
   One way to do it is to take one of the `.gt` files and do `cut -d' ' -f1 $gt_file > list.txt`
- `sample_table`: This is the 1kg super/sub populations table. (not currently in use,
   but maybe someday)
- `segments_dir`: Directory containing the `.gt` files
- `input_files`: a yaml list of the `.pos` files we are useing for the test set.
- `outdir`: output directory of all the results
- `batch_size`: batch size for model during inference.  128 is a good default.
- `num_workers`: number of processes that the model's dataloader uses during inference.
   4 is probably a good default.
- `num_neighbors`: the number of neighbors to use for k-nn queries.  I've been using 10, but
  you can experiment with other values for different analyses.

# Running without snakemake
## Setup Environment
In this directory there are two environment yaml files, depending on your hardware
setup: `torch-cpu.yaml` and `torch-gpu.yaml`.  We can create an environment by running
```
conda env create -f $YAML_FILE
```

## Running Inference
To apply the model in inference mode we run the `python encode_samples.py`
```
usage: encode_samples.py [-h] --encoder ENCODER --outdir OUTDIR --files FILES [FILES ...]
                         [--batch-size BATCH_SIZE] [--gpu] [--num-workers NUM_WORKERS]

options:
  -h, --help            show this help message and exit
  --encoder ENCODER     Path to the encoder checkpoint
  --outdir OUTDIR       Output directory (made if non-existent)
  --files FILES [FILES ...]
                        Files to encode
  --batch-size BATCH_SIZE
                        Batch size for inference
  --gpu                 Use GPU
  --num-workers NUM_WORKERS
                        Number of workers for data loading.
```

The model runs directly on positional encodings files.  Here's an example run

```
$ python encode_samples.py \
    --encoder 'swa.checkpoints/siamese-epoch=10-val_loss=0.01.ckpt' \
    --outdir '$PATH_TO_OUTDIR' \
    --files '$DATA_DIR/1KG.data.seg.9*.encoded' \
    --batch-size 128 \ # sensible default, but could be larger on bigger GPUs or CPU
    --gpu \ # set this flag if using a gpu
    --num-workers $N_PROCESSES
```
In the files list I used a wildcard to catch any files that matched that pattern,
but you could also just provide a space separated list of files there.


The output will look something like this
```
HG00096_0	90	-0.24322154 -1.4407103 -0.4442817   ...
HG00096_1	90	-0.1990201 -0.9755058 -0.43592215   ...
HG00097_0	90	-0.4929338 -1.8236687 -0.68836486   ...
HG00097_1	90	-0.3018111 -1.7684044 -0.52830446   ...
HG00099_0	90	-0.45660067 -1.7740071 -0.61605495  ...
HG00099_1	90	-0.378156 -1.8069385 -0.4742077     ...
HG00100_0	90	-0.2851319 -1.6945909 -0.5311256    ...
HG00100_1	90	-0.2686975 -1.0184398 -0.5486269    ...
HG00101_0	90	-0.20446123 -0.36143845 -0.32215124 ...
HG00101_1	90	-0.079634786 -0.7559077 -0.24679162 ...
...
```
Columns are speparated by *tab*. Column 1 is the sample name, column 2 is the
segment number and the 3rd column is a *space* separated list of floats (ie the embedding).

At the moment, embedding dimension is 512, but we can experiment with that later.

## Train a new model
- TODO
