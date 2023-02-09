
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
