#!/bin/sh

python_dir="/home/sdp/precision-medicine/python/scripts/plotting"
data_dir="/home/sdp/precision-medicine/python/scripts/plotting/txt"
png_dir="/home/sdp/precision-medicine/python/scripts/plotting/png"
avgTPenc=$datadir/avg.tp.encoding.txt
avgTPemb=$datadir/avg.tp.embedding.txt
avgPenc=$datadir/avg.p.encoding.txt
avgPemb=$datadir/avg.p.embedding.txt
medTPenc=$datadir/med.tp.encoding.txt
medTPemb=$datadir/med.tp.embedding.txt
medPenc=$datadir/med.p.encoding.txt
medPemb=$datadir/med.p.embedding.txt

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
