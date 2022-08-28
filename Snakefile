from types import SimpleNamespace
configfile: "config.yaml" # path to the config
config = SimpleNamespace(**config)

LD_LIBRARY_PATH = f"{config.conda_dir}/lib"
shell.prefix("""
set -euo pipefail;
export LD_LIBRARY_PATH=\"{LD_LIBRARY_PATH}\";
""".format(LD_LIBRARY_PATH=LD_LIBRARY_PATH))

rule all:
	input:
		f"{config.samples_dir}/sampleIDs.ALL",
		f"{config.bin_dir}/slice-vcf",
		f"{config.bin_dir}/encode-vcf", 
		f"{config.bin_dir}/faiss",
                f"{config.segments_out_dir}/segments.vcf.done", 
		f"{config.segments_out_dir}/segments.plink.done",
		f"{config.segments_out_dir}/segments.encoding.done",
		f"{config.segments_out_dir}/segments.database-embeddings.done",
		f"{config.segments_out_dir}/segments.query-embeddings.done",
		f"{config.segments_out_dir}/segments.enc_faissL2.done",
		f"{config.segments_out_dir}/segments.emb_faissL2.done"

rule sample_IDs:
	input:
		vcf=f"{config.vcf_file}"
	output:
		sample_IDs=f"{config.samples_dir}/sampleIDs.ALL"
	message:
		"Getting list of all sample IDs from VCF file..."
	shell:
		"bcftools query -l {input.vcf} > {output.sample_IDs}"


rule split_vcf_COMPILE:
	input:
        	slice_vcf=f"{config.src_dir}/slice_vcf.cpp",
		read_config=f"{config.src_dir}/read_config.cpp"
	output:
                bin=f"{config.bin_dir}/slice-vcf"
	message:
		"Compiling--slice vcf into segments..."
	shell:
                "g++ " \
		" {input.slice_vcf} " \
                " {input.read_config} " \
                " -I {config.include_dir}/" \
                " -lhts" \
                " -o {output.bin}"

rule split_vcf_EXECUTE:
	input:
		bin=f"{config.bin_dir}/slice-vcf",
		config_file=f"{config.configs_dir}/sample.config"
	output:
		done=f"{config.segments_out_dir}/segments.vcf.done"
	message: 
		"Executing--slice vcf into segments..."
	shell:
		"./{input.bin} {input.config_file} > {output.done}" \
		" && touch {output.done}"


rule plink_genome_IBD:
	input:
		vcf_done=f"{config.segments_out_dir}/segments.vcf.done"
	output:
		done=f"{config.segments_out_dir}/segments.plink.done"
	message:
		"Running plink on all encodings"
	shell:
		"for vcf_f in {config.segments_out_dir}/*.vcf; do" \
		"       filename=$(basename $vcf_f);" \
		"	seg_name=${{filename%.*}};" \
		"	plink --vcf $vcf_f --genome --out {config.segments_out_dir}/${{seg_name}} > {output.done};" \
		"done" \
		" && rm {config.segments_out_dir}/*.log" \
		" && rm {config.segments_out_dir}/*.nosex" \
		" && touch {output.done}"


rule encode_vcf_COMPILE:
	input:
		encode_vcf=f"{config.src_dir}/encode_vcf.cpp",
		read_config=f"{config.src_dir}/read_config.cpp", 
		map_encodings=f"{config.src_dir}/map_encodings.cpp",
		utils=f"{config.src_dir}/utils.cpp"
	output:
		bin=f"{config.bin_dir}/encode-vcf"
	message:
		"Compiling--encode vcf slices..."
	shell:
		"g++ " \
                " {input.encode_vcf} " \
                " {input.read_config} " \
		" {input.map_encodings} "\
		" {input.utils} " \
                " -I {config.include_dir}/" \
                " -lhts" \
                " -o {output.bin}"

rule encode_vcf_EXECUTE:
	input:
		bin=f"{config.bin_dir}/encode-vcf",
		sample_IDs=f"{config.samples_dir}/sampleIDs.ALL", 
		config_file=f"{config.configs_dir}/sample.config"
	output:
		done=f"{config.segments_out_dir}/segments.encoding.done"
	message:
		"Executing--encode vcf slices..."
	shell:
		"for vcf_f in {config.segments_out_dir}/*.vcf; do" \
		"	filename=$(basename $vcf_f);" \
                "       seg_name=${{filename%.*}};" \
		"	./{input.bin} {input.config_file} {input.sample_IDs} $vcf_f {config.segments_out_dir}/${{seg_name}}.encoding > {output.done};" \
		"done" \
		" && touch {output.done}"


rule run_model_database_data:
	input:
		encodings=f"{config.segments_out_dir}/segments.encoding.done",
		script=f"{config.model_dir}/create_vectors.py", 
		database_samples=f"{config.samples_dir}/train_samples.txt", 
		model=f"{config.model_dir}/base_model.h5",
		sample_IDs=f"{config.samples_dir}/sampleIDs.ALL",
		
	output:
		done=f"{config.segments_out_dir}/segments.database-embeddings.done"
	message:
		"Running model for all vcf slices on database data. (Getting embeddings...)"
	shell:
		"for encoding_f in {config.segments_out_dir}/*.encoding; do" \
		"       filename=$(basename $encoding_f);" \
		"       seg_name=${{filename%.*}};" \
		"	python {input.script}" \
    		" 	 --sample-list {input.database_samples}" \
    		"	 --model-path {input.model}" \
    		"	 --output-path {config.segments_out_dir}/${{seg_name}}.database-embedding" \
    		"	 --batch-size {config.batch_size}" \
    		"	 --sample_id_filename {input.sample_IDs}" \
    		"	 --genotype_filename $encoding_f;"
		"done" \
		" && touch {output.done}"

rule run_model_query_data:
	input:
		encodings=f"{config.segments_out_dir}/segments.encoding.done",
		script=f"{config.model_dir}/create_vectors.py", 
		query_samples=f"{config.samples_dir}/test_samples.txt", 
		model=f"{config.model_dir}/base_model.h5",
		sample_IDs=f"{config.samples_dir}/sampleIDs.ALL",
		
	output:
		done=f"{config.segments_out_dir}/segments.query-embeddings.done"
	message:
		"Running model for all vcf slices on query data. (Getting embeddings...)"
	shell:
		"for encoding_f in {config.segments_out_dir}/*.encoding; do" \
		"       filename=$(basename $encoding_f);" \
		"       seg_name=${{filename%.*}};" \
		"	python {input.script}" \
    		" 	 --sample-list {input.query_samples}" \
    		"	 --model-path {input.model}" \
    		"	 --output-path {config.segments_out_dir}/${{seg_name}}.query-embedding" \
    		"	 --batch-size {config.batch_size}" \
    		"	 --sample_id_filename {input.sample_IDs}" \
    		"	 --genotype_filename $encoding_f;"
		"done" \
		" && touch {output.done}"


rule faiss_COMPILE:
	input:
		f"{config.segments_out_dir}/segments.encoding.done",
                f"{config.segments_out_dir}/segments.database-embeddings.done",
                f"{config.segments_out_dir}/segments.query-embeddings.done",
		main=f"{config.src_dir}/single_faiss.cpp",
		utils=f"{config.src_dir}/utils.cpp",
		build=f"{config.src_dir}/build_index.cpp",
		search=f"{config.src_dir}/search_index.cpp",
		read=f"{config.src_dir}/read_encodings.cpp"
	output:
		bin=f"{config.bin_dir}/faiss"
	message:
		"Compiling--FAISS..."
	shell:
		"g++ " \
                " {input.main} " \
                " {input.utils} " \
                " {input.build} " \
                " {input.search} "\
		" {input.read} "\
                " -I {config.include_dir}/" \
		" -I {config.conda_dir}/include/" \
		" -L {config.conda_dir}/lib/" \
                " -lfaiss" \
                " -o {output.bin}" \

rule faiss_encoding_EXECUTE:
	input:
		encodings=f"{config.segments_out_dir}/segments.encoding.done",
		bin=f"{config.bin_dir}/faiss",
		database_IDs=f"{config.samples_dir}/train_samples.txt",
		query_IDs=f"{config.samples_dir}/test_samples.txt"
	output:
		done=f"{config.segments_out_dir}/segments.enc_faissL2.done"
	message:
		"Executing--run FAISS L2 on all input encoding segments"
	shell:
		"for f in {config.segments_out_dir}/*.encoding; do" \
		"       filename=$(basename $f);" \
		"	seg_name=${{filename%.*}};" \
		"	./{input.bin} {input.database_IDs} {config.segments_out_dir}/${{seg_name}}.encoding {input.query_IDs} {config.segments_out_dir}/${{seg_name}}.encoding {config.k} {config.enc_delim} > {config.segments_out_dir}/${{seg_name}}.enc_faissL2;" \ 
		"done" \
		" && touch {output.done}"


rule faiss_embedding_EXECUTE:
	input:
		db_embeddings=f"{config.segments_out_dir}/segments.database-embeddings.done",
		q_embeddings=f"{config.segments_out_dir}/segments.query-embeddings.done",
		bin=f"{config.bin_dir}/faiss",
		database_IDs=f"{config.samples_dir}/train_samples.txt",
		query_IDs=f"{config.samples_dir}/test_samples.txt"
	output:
		done=f"{config.segments_out_dir}/segments.emb_faissL2.done"
	message:
		"Executing--run FAISS L2 on all input embedding segments"
	shell:
		"for f in {config.segments_out_dir}/*.encoding; do" \
		"       filename=$(basename $f);" \
		"	seg_name=${{filename%.*}};" \
		"	./{input.bin} {input.database_IDs} {config.segments_out_dir}/${{seg_name}}.database-embedding {input.query_IDs} {config.segments_out_dir}/${{seg_name}}.query-embedding {config.k} {config.emb_delim} > {config.segments_out_dir}/${{seg_name}}.emb_faissL2;" \ 
		"done" \
		" && touch {output.done}"
