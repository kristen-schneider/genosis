#!/bin/env bash
# evaluate test set with a gpu node on fiji
# SBATCH -p titan # TODO check if this is correct
# SBATCH --job-name=eval_test_set
# SBATCH --gres=gpu:0
# SBATCH --ntasks=1
# SBATCH --nodes=1
# SBATCH --time=4:00:00
# SBATCH --mem=16gb
# SBATCH --output="./eval_test_set.out"
# SBATCH --error="./eval_test_set.err"
# --------------------------------------------------------------

data_dir=/scratch/Shares/layer/projects/toy_model_data
python evaluation_loop.py \
    --model_path $data_dir/base_model.h5 \
    --test_samples $data_dir/testing.txt \
    --sample_ids $data_dir/ALL.sampleIDs \
    --sample_pairs $data_dir/chr0.seg.0.cnn \
    --genotypes $data_dir/chr0.seg.0.encoding \
    --batch_size 64

