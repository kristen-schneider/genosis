from types import SimpleNamespace
configfile: "config.yaml" # path to the config
config = SimpleNamespace(**config)

rule all:
	input:
		f"{config.segments_out_dir}/segments.encoding.done",		
		f"{config.segments_out_dir}/segments.plink.done",		
		f"{config.segments_out_dir}/segments.faissL2.done",
		f"{config.segments_out_dir}/segments.cnn.done"	

rule split_encode_vcf_COMPILE:
	output:
		bin=f"{config.bin_dir}/segment"
	message:
		"Compiling--slice vcf into segments and encode"
	shell:
		"g++ {config.src_dir}/main.cpp" \
		" {config.src_dir}/read_config.cpp" \
		" {config.src_dir}/map_encodings.cpp" \
		" {config.src_dir}/slice_vcf.cpp" \
		" {config.src_dir}/encode_vcf.cpp" \
		" -I {config.include_dir}/" \
		" -lhts" \
		" -o {output.bin}"

rule split_encode_vcf_EXECUTE:
	input:
		bin=f"{config.bin_dir}/segment",
		config_file=f"{config.configs_dir}/segment_config"
	output:
		done=f"{config.segments_out_dir}/segments.encoding.done"
	message: 
		"Executing--slice vcf into segments and encode"
	shell:
		"./{input.bin} {input.config_file}" \
		" && touch {output.done}"

rule plink_genome_IBD:
	input:
		encoding_done=f"{config.segments_out_dir}/segments.encoding.done"
	output:
		done=f"{config.segments_out_dir}/segments.plink.done"
	message:
		"Running plink on all encodings"
	shell:
		"for vcf_f in {config.segments_out_dir}/*.vcf; do" \
		"	plink --vcf $vcf_f --genome;" \
		"       filename=$(basename $vcf_f);" \
		"	seg_name=${{filename%.*}};" \
		"	mv plink.genome {config.segments_out_dir}/${{seg_name}}.genome;" \
		"done" \
		" && touch {output.done}"

rule faiss_L2_COMPILE:
	input:
		encodings=f"{config.segments_out_dir}/segments.encoding.done"
	output:
		bin=f"{config.bin_dir}/single_faiss"
	message:
		"Compiling--run FAISS L2 on all input segments"
	shell:
		"g++ {config.src_dir}/single_faiss.cpp" \
		" {config.src_dir}/buildIndex.cpp" \
		" {config.src_dir}/readEncoding.cpp" \
		" {config.src_dir}/searchIndex.cpp" \
		" {config.src_dir}/utils.cpp" \
		" -I {config.include_dir}/" \
		" -I {config.conda_dir}/include/" \
		" -L {config.conda_dir}/lib/" \
		" -lfaiss" \
		" -o {output.bin}"		

rule faiss_L2_EXECUTE:
	input:
		bin=f"{config.bin_dir}/single_faiss",
	output:
		done=f"{config.segments_out_dir}/segments.faissL2.done"
	message:
		"Executing--run FAISS L2 on all input segments"
	shell:
		"for encoding_f in {config.segments_out_dir}/*.encoding; do" \
		"       filename=$(basename $encoding_f);" \
		"	seg_name=${{filename%.*}};" \
		"	./{input.bin} $encoding_f $encoding_f {config.k} {config.delim}> {config.segments_out_dir}/${{seg_name}}.faissL2;" \ 
		"done" \
		" && touch {output.done}"


rule generate_cnn_input_file:
	input:
		sample_IDs=f"{config.sample_IDs}",
	output:
		done=f"{config.segments_out_dir}/segments.cnn.done"
	message:
		"Generating CNN input files for all segments"
	shell:
		"for encoding_f in {config.segments_out_dir}/*.encoding; do" \
		"	filename=$(basename $encoding_f);" \
                "       seg_name=${{filename%.*}};" \
		"	python {config.python_dir}/scripts/cnn/make_CNN_input.py {input.sample_IDs} $encoding_f {config.segments_out_dir}/${{seg_name}}.genome {config.segments_out_dir}/${{seg_name}}.cnn;" \
		"done" \
		" && touch {output.done}"
