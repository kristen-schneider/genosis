example_dir="/Users/krsc0813/precision-medicine/example/"
cd $example_dir

rm sample_IDs.txt
rm vcf_bp.txt
rm segment_boundary.map
rm interpolated.map
rm *hap_IDs.txt

rm err/*
rm out/*
rm -r log/
rm -r vcf_segments/
rm -r encodings/
rm -r embeddings/
rm -r svs_index/
rm -r svs_results/
rm TOP_HITS.txt
#rm -r svs_sample_results/
