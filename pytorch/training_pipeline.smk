import os
import random
import sys
from glob import glob
from types import SimpleNamespace

# ==============================================================================
# Project Configuration
# ==============================================================================
configfile: "conf/exp_config.yaml"
config = SimpleNamespace(**config)


segments = [
  x.split(config.file_segment_delimiter)[config.seg_offset]
  for x in glob(f'{config.segments_dir}/*.{config.gt_ext}')
  if x.split(config.file_segment_delimiter)[config.seg_offset] != '168' # lol
]
if segments == []:
  raise ValueError(
    f"No segments found\n"
    f"{config.segments_dir=}\n"
    f"{config.gt_ext=}\n"
    f"{config.file_segment_delimiter=}\n"
    f"{config.seg_offset=}"
  )

# ==============================================================================
# train, val, test split by segment
# ==============================================================================
random.seed(config.random_seed)
shuffled_segments = random.sample(segments, len(segments))

train_ratio = 0.8 # TODO make this a config option
val_ratio = test_ratio = (1.0 - train_ratio) / 2.0

train_segments = shuffled_segments[
  : int(train_ratio * len(segments))
]
val_segments = shuffled_segments[
  int(train_ratio * len(segments)) : int((train_ratio + val_ratio) * len(segments))
]
test_segments = shuffled_segments[
  int((train_ratio + val_ratio) * len(segments)) :
]

# any empty sets?
if not train_segments or not val_segments or not test_segments:
  raise ValueError(
    f"Train, val, test split failed\n"
    f"{train_ratio=}\n"
    f"{val_ratio=}\n"
    f"{test_ratio=}\n"
    f"{train_segments=}\n"
    f"{val_segments=}\n"
    f"{test_segments=}\n"
  )

# check if the sets are disjoint
if len(set(train_segments) & set(val_segments) & set(test_segments)) > 0:
  raise ValueError(
    f"Train/val/test sets are not disjoint\n"
    f"{train_segments=}\n"
    f"{val_segments=}\n"
    f"{test_segments=}\n"
  )


# ==============================================================================
# Rules
# ==============================================================================
rule All:
  input:
    # trained model
    directory(f"{config.outdir}/{config.model_prefix}.checkpoints"),

    # train memmaps
    f"{config.outdir}/training_set/P1.mmap",
    f"{config.outdir}/training_set/P2.mmap",

    # val memmaps
    f"{config.outdir}/validation_set/P1.mmap",
    f"{config.outdir}/validation_set/P2.mmap",

    # train/test/split segement numbers
    f"{config.outdir}/train.segments",
    f"{config.outdir}/val.segments",
    f"{config.outdir}/test.segments",

    # training set
    expand(f"{config.outdir}/training_set/D.txt", segment=train_segments),

    # validation set
    expand(f"{config.outdir}/validation_set/D.txt", segment=val_segments),

    # stats
    f"{config.outdir}/bin_sampled_distribution.pdf",
    f"{config.outdir}/distribution.pdf",

rule SampleSubpops:
  input:
    config.sample_table,
  output:
    f"{config.outdir}/subpops/pairings.{{segment}}.txt"
  shell:
    f"""
    python exploration/sample_subpops.py \
      --sample_table {{input}} \
      --samples_list {config.samples_list} \
      --output {{output}} \
      --N {config.num_pairings}
    """

rule MakeGTMemmap:
  """
  For each segment, make a numpy memmap of the genotype matrix
  """
  input:
    f"{config.segments_dir}/{config.segment_prefix}.{{segment}}.{config.gt_ext}",
  output:
    f"{config.outdir}/gt_mmap/segment.{{segment}}.mmap"
  conda:
    "envs/mmap_ninja.yaml"
  shell:
    f"""
    python exploration/make_gt_mmap.py \
      --gts {{input}} \
      --output {{output}} \
    """

rule ComputeDistances:
  """
  For each sample, get the distances of all the samples paired with it.
  Since the paired samples don't have haplotype numbers, randomly choose one
  """
  input:
    memmap = rules.MakeGTMemmap.output,
    pairings = rules.SampleSubpops.output,
  output:
    f"{config.outdir}/distances/distances.{{segment}}.txt"
  shell:
    f"""
    python exploration/compute_distances.py \
      --memmap {{input.memmap}} \
      --pairings {{input.pairings}} \
      --output {{output}} \
      --samples_list {config.samples_list}
    """

rule PlotDistribution:
  """
  Plot the distribution of distances
  """
  input:
    rules.ComputeDistances.output
  output:
    temp(f"{config.outdir}/plots/distribution.{{segment}}.png")
  shell:
    f"""
    python exploration/plot_distribution.py \
      --segment {{wildcards.segment}} \
      --distances {{input}} \
      --output {{output}}
    """

rule CatImagesPDF:
  """
  Concatenate all the images into a single pdf
  """
  input:
    expand(rules.PlotDistribution.output, segment=segments)
  output:
    f"{config.outdir}/distribution.pdf"
  shell:
    f"""
    convert {{input}} {{output}}
    """

rule BinSampling:
  """
  For each segment, divide the space of possible distances into bins and
  sample from each bin to get a representative sample of distances.
  """
  input:
    memmap = rules.MakeGTMemmap.output,
    distances = rules.ComputeDistances.output,
  output:
    f"{config.outdir}/bin_sampling/segment.{{segment}}.txt"
  shell:
    f"""
    python exploration/bin_sampling.py \
      --distances {{input.distances}} \
      --output {{output}} \
      --num_bins 20 \
      --max_samples 500
    """

rule PlotResampledDistribution:
  """
  Plot the distribution of distances from the resampled distances
  """
  input:
    rules.BinSampling.output
  output:
    temp(f"{config.outdir}/resampled_plots/distribution.{{segment}}.png")
  shell:
    f"""
    python exploration/plot_distribution.py \
      --segment {{wildcards.segment}} \
      --distances {{input}} \
      --output {{output}}
    """

rule CatResampledImagesPDF:
  """
  Concatenate all the images into a single pdf
  """
  input:
    expand(rules.PlotResampledDistribution.output, segment=segments)
  output:
    f"{config.outdir}/bin_sampled_distribution.pdf"
  shell:
    f"""
    convert {{input}} {{output}}
    """

rule WriteTrainTestSplit:
  """
  Write the train and test segment numbers to a file
  """
  output:
    train = f"{config.outdir}/train.segments",
    val   = f"{config.outdir}/val.segments",
    test  = f"{config.outdir}/test.segments"
  run:
    with open(output.train, 'w') as f:
      for segment in sorted(map(int, train_segments)):
        f.write(f"{segment}\n")
    with open(output.val, 'w') as f:
      for segment in sorted(map(int, val_segments)):
        f.write(f"{segment}\n")
    with open(output.test, 'w') as f:
      for segment in sorted(map(int, test_segments)):
        f.write(f"{segment}\n")

rule MakeTrainingSet:
  """
  From the resampled distances, make a training set of pairs
  of samples, over all segments.
  """
  input:
    pos_files = expand(
      f"{config.segments_dir}/{config.segment_prefix}.{{segment}}.{config.pos_ext}",
      segment=train_segments
    ),
    distances = expand(rules.BinSampling.output, segment=train_segments),
  output:
    P1 = temp(f"{config.outdir}/training_set/P1.txt"),
    P2 = temp(f"{config.outdir}/training_set/P2.txt"),
    D = f"{config.outdir}/training_set/D.txt",

  shell:
    f"""
    python exploration/make_dataset.py \
      --pos_files {{input.pos_files}} \
      --distance_files {{input.distances}} \
      --P1 {{output.P1}} \
      --P2 {{output.P2}} \
      --D {{output.D}} \
    """

rule MakeValidationSet:
  """
  From the resampled distances, make a validation set of pairs
  of samples, over all segments.
  """
  input:
    pos_files = expand(
      f"{config.segments_dir}/{config.segment_prefix}.{{segment}}.{config.pos_ext}",
      segment=val_segments
    ),
    distances = expand(rules.BinSampling.output, segment=val_segments),
  output:
    P1 = temp(f"{config.outdir}/validation_set/P1.txt"),
    P2 = temp(f"{config.outdir}/validation_set/P2.txt"),
    D = f"{config.outdir}/validation_set/D.txt",

  shell:
    f"""
    python exploration/make_dataset.py \
      --pos_files {{input.pos_files}} \
      --distance_files {{input.distances}} \
      --P1 {{output.P1}} \
      --P2 {{output.P2}} \
      --D {{output.D}} \
    """

rule MakeTrainMmaps:
  """
  Make memmaps for the training set position vectors
  """
  input:
    P1 = rules.MakeTrainingSet.output.P1,
    P2 = rules.MakeTrainingSet.output.P2,
  output:
    P1 = directory(f"{config.outdir}/training_set/P1.mmap"),
    P2 = directory(f"{config.outdir}/training_set/P2.mmap")
  conda:
    "envs/mmap_ninja.yaml"
  shell:
    f"""
    python exploration/generate_mmaps.py \
      --inP1 {{input.P1}} \
      --inP2 {{input.P2}} \
      --outP1 {{output.P1}} \
      --outP2 {{output.P2}}
    """

rule MakeValMmaps:
  """
  Make memmaps for the validation set position vectors
  """
  input:
    P1 = rules.MakeValidationSet.output.P1,
    P2 = rules.MakeValidationSet.output.P2,
  output:
    P1 = directory(f"{config.outdir}/validation_set/P1.mmap"),
    P2 = directory(f"{config.outdir}/validation_set/P2.mmap")
  conda:
    "envs/mmap_ninja.yaml"
  shell:
    f"""
    python exploration/generate_mmaps.py \
      --inP1 {{input.P1}} \
      --inP2 {{input.P2}} \
      --outP1 {{output.P1}} \
      --outP2 {{output.P2}}
    """

rule TrainModel:
  input:
    train_segments = f"{config.outdir}/train.segments",
    val_segments = f"{config.outdir}/val.segments",
    test_segments = f"{config.outdir}/test.segments",
    P1_train = rules.MakeTrainMmaps.output.P1,
    P2_train = rules.MakeTrainMmaps.output.P2,
    P1_val = rules.MakeValMmaps.output.P1,
    P2_val = rules.MakeValMmaps.output.P2,
    D_train = rules.MakeTrainingSet.output.D,
    D_val = rules.MakeValidationSet.output.D
  output:
    model_checkpoints = directory(f"{config.outdir}/{config.model_prefix}.checkpoints")
  threads:
    config.n_workers
  conda:
    'envs/torch-gpu.yaml'
  shell:
    f"""
    python train_model.py \
      --outdir {config.outdir} \
      --model_prefix {config.model_prefix} \
      --P1_train {{input.P1_train}} \
      --P2_train {{input.P2_train}} \
      --P1_val {{input.P1_val}} \
      --P2_val {{input.P2_val}} \
      --D_train {{input.D_train}} \
      --D_val {{input.D_val}} \
      --model_config {config.model_config}
    """


