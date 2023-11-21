#!/usr/bin/env bash
#SBATCH --partition=short
#SBATCH --job-name=plot
#SBATCH --output=./out/plot.out
#SBATCH --error=./err/plot.err
#SBATCH --time=0-23:00:00
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mail-type=NONE
#SBATCH --mail-user=krsc0813@colorado.edu

biobagss_dir="/Users/krsc0813/biobagg_analysis/"
pop_file="/Users/krsc0813/1KG_data/samples.ancestry"
trio_file="/Users/krsc0813/1KG_data/samples.trios"
out_dir="/scratch/Users/krsc0813/chr1_22/"
knn_file=$out_dir"TOP_HITS.txt"

python $biobagss_dir"plotting/evaluate_ancestry" \
  --pop $pop_file\
  --knn $knn_file\
  --png $out_dir"violin_plots/"

python $biobagss_dir"plotting/search_time.py" \
  --num_samples 3202 \
  --in_file $out_dir"log/search.log" \
  --out_file $out_dir"time_plots/1kg_search_time_histo.png" \
  --height 4 \
  --width 8

python $biobagss_dir"plotting/top_hits_umap.py" \
  --in_file $knn_file \
  --label_file $pop_file \
  --out_file $out_dir"other_plots/top_hit_umap.png"

python plotting/top_hits_pca.py \
  --in_file $knn_file \
  --label_file $pop_file \
  --out_file $out_dir"other_plots/top_hit_pca.png"

python plotting/plot_topk_score_distro.py \
  --topk_file $knn_file \
  --label_file $pop_file \
  --out_file $out_dir"other_plots/topk_score_distro.png"  \
  --height 2 \
  --width 8 \
  --bins 100 \
  --ped_file $trio_file

python plotting/plot_topk_score_distro.py \
  --topk_file $knn_file \
  --label_file $pop_file \
  --out_file $out_dir"other_plots/topk_score_distro.log.png"  \
  --height 2 \
  --width 8 \
  --bins 100 \
  --ped_file $trio_file \
  --log


