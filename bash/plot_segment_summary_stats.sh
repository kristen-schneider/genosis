#!/bin/sh

python_dir="/home/sdp/precision-medicine/python/scripts/plotting"
data_dir="/home/sdp/precision-medicine/python/scripts/plotting/txt"
png_dir="/home/sdp/precision-medicine/python/scripts/plotting/png"
avgTPenc=$data_dir"/avg.tp.encoding.txt"
avgTPemb=$data_dir"/avg.tp.embedding.txt"
avgPenc=$data_dir"/avg.p.encoding.txt"
avgPemb=$data_dir"/avg.p.embedding.txt"
medTPenc=$data_dir"/med.tp.encoding.txt"
medTPemb=$data_dir"/med.tp.embedding.txt"
medPenc=$data_dir"/med.p.encoding.txt"
medPemb=$data_dir"/med.p.embedding.txt"

python $python_dir/plot_summary_segments.py \
 --avg_tp_encoding $avgTPenc \
 --avg_tp_embedding $avgTPemb \
 --avg_p_encoding $avgPenc \
 --avg_p_embedding $avgPemb \
 --med_tp_encoding $medTPenc \
 --med_tp_embedding $medTPemb \
 --med_p_encoding $medPenc \
 --med_p_embedding $medPemb \
 --avg_png $png_dir/avg.png \
 --med_png $png_dir/med.png
