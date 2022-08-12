#configfile: "config.yaml"

rule all:
	input:
		"setup.done",
		"data/segments/seg_5000/segments.encoding.done", 
		"data/segments/seg_5000/segments.genome.done",
		"data/segments/seg_5000/segments.faissL2.done"

rule set_up:
	input:
		condalib="/home/sdp/miniconda3/envs/precision-medicine/lib"
	output:
		done="setup.done"
	message: 
		"Setting up environment"
	shell:
		" source ~/miniconda3/etc/profile.d/conda.sh" \
		" && conda activate precision-medicine" \
		" && PATH=$PATH:~/plink_linux_x86_64_20220402/" \
		" && export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:{input.condalib}" \
		" && touch {output.done}"
		

rule split_encode_vcf_COMPILE:
	input:
		main="cpp/src/main.cpp",
		read_config="cpp/src/read_config.cpp",
		map_encodings="cpp/src/map_encodings.cpp",
		slice_vcf="cpp/src/slice_vcf.cpp",
		encode_vcf="cpp/src/encode_vcf.cpp",
		include="cpp/include/"
	output:
		bin="cpp/bin/segment"
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
		bin="cpp/bin/segment",
		config_file="cpp/configs/segment_config"
	output:
		done="data/segments/seg_5000/segments.encoding.done"
	message: 
		"Executing--slice vcf into segments and encode"
	shell:
		"./{input.bin} {input.config_file}" \
		" && touch {output.done}"


rule plink_genome_IBD:
	input:
		encoding_done="data/segments/seg_5000/segments.encoding.done"
	output:
		done="data/segments/seg_5000/segments.genome.done"
	message:
		"Running plink on all encodings"
	shell:
		"bash bash/run_plink.sh" \
		" && touch {output.done}"


rule faiss_L2_COMPILE:
	input:
		main="cpp/src/single_faiss.cpp",
		build="cpp/src/buildIndex.cpp",
		read="cpp/src/readEncoding.cpp",
		search="cpp/src/searchIndex.cpp",
		utils="cpp/src/utils.cpp",
		include="cpp/include/",
		condainclude="/home/sdp/miniconda3/envs/precision-medicine/include",
		condalib="/home/sdp/miniconda3/envs/precision-medicine/lib"
	output:
		bin="cpp/bin/single_faiss"
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
		bin="cpp/bin/single_faiss",
		data_dir="data/segments/seg_5000/"
	output:
		done="data/segments/seg_5000/segments.faissL2.done"
	message:
		"Executing--run FAISS L2 on all input segments"
	shell:
		"for encoding_f in {input.data_dir}*encoding; do" \
		"       filename=$(basename $encoding_f);" \
		"	seg_name=${{filename%.*}};"
		#"	echo ./{input.bin} $encoding_f $encoding_f;"
		#"	./{input.bin} $encoding_f $encoding_f 2548;" \ 
		"	./{input.bin} $encoding_f $encoding_f 2548 > {input.data_dir}$seg_name'.faissL2';" \ 
		"done" \
		" && touch {output.done}"
