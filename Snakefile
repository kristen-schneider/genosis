from types import SimpleNamespace
configfile: "config.yaml" # path to the config
config = SimpleNamespace(**config)

rule all:
	input:
		f"{config.segments_out_dir}/segments.encoding.done",		
		f"{config.segments_out_dir}/segments.genome.done",		
		f"{config.segments_out_dir}/segments.faissL2.done",
		f"{config.segments_out_dir}/segments.cnn.done"	

rule split_encode_vcf_COMPILE:
	input:
		setup="setup.done",
		main=f"{config.src_dir}/main.cpp",
		read_config=f"{config.src_dir}/read_config.cpp",
		map_encodings=f"{config.src_dir}/map_encodings.cpp",
		slice_vcf=f"{config.src_dir}/slice_vcf.cpp",
		encode_vcf=f"{config.src_dir}/encode_vcf.cpp",
		include=f"{config.include_dir}/"
	output:
		bin=f"{config.bin_dir}/segment"
	message:
		"Compiling--slice vcf into segments and encode"
	shell:
		"g++ {input.main}" \
			" {input.read_config}" \
			" {input.map_encodings}" \
			" {input.slice_vcf}" \
			" {input.encode_vcf}" \
			" -I {input.include}" \
			" -lhts" \
			" -o {output.bin}"
		

rule split_encode_vcf_EXECUTE:
	input:
		setup="setup.done",
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
		setup="setup.done",
		data_dir=f"{config.segments_out_dir}/",
		encoding_done=f"{config.segments_out_dir}/segments.encoding.done"
	output:
		done=f"{config.segments_out_dir}/segments.genome.done"
	message:
		"Running plink on all encodings"
	shell:
		"for vcf_f in {input.data_dir}*.vcf; do" \
		"	echo $vcf_f;" \
		"	plink --vcf $vcf_f --genome;" \
		"       filename=$(basename $vcf_f);" \
		"	seg_name=${{filename%.*}};" \
		"	mv plink.genome {input.data_dir}${{seg_name}}.genome;" \
		"done" \
		" && touch {output.done}"
		

rule faiss_L2_COMPILE:
	input:
		setup="setup.done",
		encoding_done=f"{config.segments_out_dir}/segments.encoding.done",
		main=f"{config.src_dir}/single_faiss.cpp",
		build=f"{config.src_dir}/buildIndex.cpp",
		read=f"{config.src_dir}/readEncoding.cpp",
		search=f"{config.src_dir}/searchIndex.cpp",
		utils=f"{config.src_dir}/utils.cpp",
		include=f"{config.include_dir}",
		condainclude=f"{config.conda_dir}/include",
		condalib=f"{config.conda_dir}/lib"
	output:
		bin=f"{config.bin_dir}/single_faiss"
	message:
		"Compiling--run FAISS L2 on all input segments"
	shell:
		"g++ {input.main}" \
                        " {input.build}" \
                        " {input.read}" \
                        " {input.search}" \
                        " {input.utils}" \
                        " -I {input.include}" \
                        " -I {input.condainclude}" \
			" -L {input.condalib}" \
			" -lfaiss" \
                        " -o {output.bin}"

rule faiss_L2_EXECUTE:
	input:
		setup="setup.done",
		bin=f"{config.bin_dir}/single_faiss",
		data_dir=f"{config.segments_out_dir}/"
	output:
		done=f"{config.segments_out_dir}/segments.faissL2.done"
	message:
		"Executing--run FAISS L2 on all input segments"
	shell:
		"for encoding_f in {input.data_dir}*.encoding; do" \
		"       filename=$(basename $encoding_f);" \
		"	seg_name=${{filename%.*}};"
		"	./{input.bin} $encoding_f $encoding_f {config.k} > {input.data_dir}$seg_name'.faissL2';" \ 
		"done" \
		" && touch {output.done}"


rule generate_cnn_input_file:
	input:
		setup="setup.done",
		sample_IDs=f"{config.sample_IDs}",
		data_dir=f"{config.segments_out_dir}/", 
		make_CNN=f"{config.python_dir}/scripts/cnn/make_CNN_input.py"
	output:
		done=f"{config.segments_out_dir}/segments.cnn.done"
	message:
		"Generating CNN input files for all segments"
	shell:
		"for encoding_f in {input.data_dir}*.encoding; do" \
		"	filename=$(basename $encoding_f);" \
                "       seg_name=${{filename%.*}};" \
		"	python {input.make_CNN} {input.sample_IDs} $encoding_f {input.data_dir}$seg_name'.genome' {input.data_dir}$seg_name'.cnn';" \
		"done" \
		" && touch {output.done}"
