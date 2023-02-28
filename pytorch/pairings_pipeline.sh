#!/bin/env bash
snakemake -s exploration_pipeline.smk -j 4 #--use-conda --conda-frontend mamba
