import os
import sys
from glob import glob
from types import SimpleNamespace

configfile: "conf/eval_config.yaml"
config = SimpleNamespace(**config)


segments = [x.split('.')[-2] for x in glob(f'{config.segments_dir}/*.encoded')][:5]

# choose N pairings for each sample from each subpopulation
# TODO put in config
N = 1

rule All:
  input:
    expand(f"{config.outdir}/subpops/pairings.{{segment}}.txt", segment=segments)

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
      --N {N}
    """


rule ComputeDistances:
  """
  For each sample, get the distances of all the samples paired with it.
  Since the paired samples don't have haplotype numbers, randomly choose one
  """
