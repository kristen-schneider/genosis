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

python $biobagss_dir"plotting/evaluate_ancestry.py" \
  --pop $pop_file\
  --knn $knn_file\
  --png $out_dir"violin_plots/"

python $biobagss_dir"plotting/search_time.py" \
  --num_samples 3202 \
  --in_file $out_dir"log/search.log" \
  --out_file $out_dir"time_plots/1kg_search_time_histo.png" \
  --height 4 \
  --width 8

#python $biobagss_dir"plotting/top_hits_umap.py" \
#  --in_file $knn_file \
#  --label_file $pop_file \
#  --out_file $out_dir"other_plots/top_hit_umap.png"

python $biobagss_dir"plotting/plot_distros.py" \
  --topk_file $knn_file \
  --pairs_file $out_dir"plink2-kin/plink2.kin0" \
  --label_file $biobagss_dir"data/igsr-1000 genomes 30x on grch38.tsv" \
  --out_file $out_dir"other_plots/distros_plink2-kin.png" \
  --width 2 \
  --height 3 \
  --bins 50 \
  --alpha 1 \
  --outlier 0.2 \
  --ped_file $trio_file \
  --x_label "Plink kinship"

python $biobagss_dir"plotting/plot_distros.py" \
  --topk_file $knn_file \
  --pairs_file $out_dir"plink-genome/plink-genome.genome" \
  --label_file $biobagss_dir"data/igsr-1000 genomes 30x on grch38.tsv" \
  --out_file $out_dir"other_plots/distros_plink-genome.png" \
  --width 2 \
  --height 3 \
  --bins 50 \
  --alpha 1 \
  --outlier 0.45 \
  --ped_file $trio_file \
  --x_label "Plink PI_HAT"

#python plotting/plot_distros.py \
#  --topk_file $knn_file \
#  --pairs_file  \
#  --label_file $biobagss_dir"data/igsr-1000 genomes 30x on grch38.tsv" \
#  --out_file $out_dir"other_plots/distros_iLASH.png" \
#  --width 2 \
#  --height 3 \
#  --bins 50 \
#  --alpha 1 \
#  --outlier 0.45 \
#  --ped_file $trio_file \
#  --x_label "iLASH IBD"
#
#python $biobagss_dir"plotting/top_hits_pca.py" \
#  --in_file $knn_file \
#  --label_file $biobagss_dir"data/igsr-1000 genomes 30x on grch38.tsv" \
#  --out_file $out_dir"other_plots/top_hit_pca.png"
#
#python $biobagss_dir"plotting/plot_topk_score_distro.py" \
#  --topk_file $knn_file \
#  --label_file $biobagss_dir"data/igsr-1000 genomes 30x on grch38.tsv" \
#  --out_file $out_dir"other_plots/topk_score_distro.png"  \
#  --height 2 \
#  --width 8 \
#  --bins 100 \
#  --ped_file $trio_file
#
#python $biobagss_dir"plotting/plot_topk_score_distro.py" \
#  --topk_file $knn_file \
#  --label_file $biobagss_dir"data/igsr-1000 genomes 30x on grch38.tsv" \
#  --out_file $out_dir"other_plots/topk_score_distro.log.png"  \
#  --height 2 \
#  --width 8 \
#  --bins 100 \
#  --ped_file $trio_file \
#  --log
#

