from types import SimpleNamespace
configfile: "/home/sdp/pmed-local/data/1KG/config_snakemake.yaml" # path to the config
#configfile: "/home/sdp/precision-medicine/example/config_snakemake.yaml" # path to the config
config = SimpleNamespace(**config)

LD_LIBRARY_PATH = f"{config.conda_dir}/lib"
shell.prefix("""
set -euo pipefail;
export LD_LIBRARY_PATH=\"{LD_LIBRARY_PATH}\";
""".format(LD_LIBRARY_PATH=LD_LIBRARY_PATH))


rule all:
	input:
		f"{config.vcf_file}",
		f"{config.data_dir}samples_IDs.txt",
		f"{config.data_dir}samples.log",
		f"{config.data_dir}plink_map.map",
		f"{config.data_dir}interpolated.map",
		f"{config.cpp_bin_dir}segment-boundary",
		f"{config.data_dir}segment_boundary.log",
                f"{config.data_dir}segment_boundary.map",
		f"{config.out_dir}slice.log",
		f"{config.cpp_bin_dir}encode",
		f"{config.out_dir}encode.log",
		#f"{config.cpp_bin_dir}faiss-l2-build",
		#f"{config.out_dir}faiss.log"
		#f"{config.data_dir}plink.log",
		#f"{config.out_dir}distance.log",
		#f"{config.data_dir}aggregate.log",
		#f"{config.data_dir}hapID.log",
		#f"{config.data_dir}all_hap_IDs.txt",	
	
# 0.1 create a file with all sample IDs
# one line per sample ID
rule get_sample_IDs:
	input:
		vcf=f"{config.vcf_file}"
	output:
		sample_IDs=f"{config.data_dir}samples_IDs.txt",
		sample_IDs_done=f"{config.data_dir}samples.log"
	message:
		"Creating a list of all sample IDs from VCF file..."
	shell:
		"bcftools query -l {input.vcf} > {output.sample_IDs}"
		" && touch {output.sample_IDs_done}"

# 0.2 create an map file with plink
rule create_map:
	input:
		vcf=f"{config.vcf_file}"
	output:
		plink_map=f"{config.data_dir}plink_map.map"
	message:
		"Using plink recode to create a map file."
	shell:
		"plink --vcf {input.vcf}" \
		" --recode 01" \
		" --output-missing-genotype ." \
		" --vcf-half-call m" \
		" --out {config.data_dir}plink_map"
	
# 0.3 create an interpolated map file
rule interpolate_map:
	input:
		plink_map=f"{config.data_dir}plink_map.map"
	output:
		interpolated_map=f"{config.data_dir}interpolated.map"
	message:
		"Using ilash_analyzer to create a new map file"
	shell:
		"python {config.python_dir}utils/interpolate_map.py" \
		" --map {input.plink_map}" \
		" --ref_map {config.ref_map}" \
		" --out_map {output.interpolated_map}"	

# 1.1 write segment boundary file (compile)
rule segment_boundary_file_compile:
	input:
		vcf_file=f"{config.vcf_file}",
		main_segment_cpp=f"{config.cpp_src_dir}main_segment.cpp",
		segment_boundary_map_cpp=f"{config.cpp_src_dir}segment_boundary_map.cpp",
		read_config_cpp=f"{config.cpp_src_dir}read_config.cpp",
		read_map_cpp=f"{config.cpp_src_dir}read_map.cpp",
		utils_cpp=f"{config.cpp_src_dir}utils.cpp"
	output:
		bin=f"{config.cpp_bin_dir}segment-boundary"
	message:
		"Compiling--write segment boundary file..."
	shell:
		"g++" \
		" {input.main_segment_cpp}" \
		" {input.segment_boundary_map_cpp}" \
		" {input.read_config_cpp}" \
		" {input.read_map_cpp}" \
		" {input.utils_cpp}" \
		" -I {config.cpp_include_dir}" \
		" -I {config.htslib_dir}" \
		" -lhts" \
		" -o {output.bin}"

# 1.2 write segment boundary file (execute)
rule segment_boundary_file_execute:
	input:
		bin=f"{config.cpp_bin_dir}segment-boundary",
		config_file=f"{config.config_file}"
	output:
		segment_boundary_log=f"{config.data_dir}segment_boundary.log",
		segment_boundary_file=f"{config.data_dir}segment_boundary.map"
	message:
		"Executing--write segment boundary file..."
	shell:
		"./{input.bin} {input.config_file} > {output.segment_boundary_log}"

# 1.3 slice VCF into segments
rule slice_VCF:
	input:
		vcf_file=f"{config.vcf_file}",
		segment_boundary_file=f"{config.data_dir}segment_boundary.map"
	output:
		slice_log=f"{config.out_dir}slice.log"
	message:
		"Slicing VCF into segments..."
	shell:
		"echo 2. ---SLICING VCF INTO SEGMENTS---;" \
		"while IFs= read -r segment start_bp end_bp; do" \
		"	echo slicing segment ${{segment}} >> {output.slice_log};" \
		"	bcftools view -h {input.vcf_file} > {config.out_dir}segment.${{segment}}.vcf;" \
		"	tabix {input.vcf_file} chr8:${{start_bp}}-${{end_bp}} >> {config.out_dir}segment.${{segment}}.vcf;" \
		" done < {input.segment_boundary_file};" \

# 2.1 encode genoypes for VCF segments (compile)
rule encode_compile:
	input:
		slice_log=f"{config.out_dir}slice.log",
		main_encode_cpp=f"{config.cpp_src_dir}main_encode.cpp",
                encode_segment_cpp=f"{config.cpp_src_dir}encode_segment.cpp",
                read_config_cpp=f"{config.cpp_src_dir}read_config.cpp",
                read_map_cpp=f"{config.cpp_src_dir}read_map.cpp",
                map_encodings_cpp=f"{config.cpp_src_dir}map_encodings.cpp",
                utils_cpp=f"{config.cpp_src_dir}utils.cpp"
	output:
		bin=f"{config.cpp_bin_dir}encode"
	message:
		"Compiling--encoding segments..."
	shell:
		"g++" \
		" {input.main_encode_cpp}" \
		" {input.encode_segment_cpp}" \
		" {input.read_config_cpp}" \
		" {input.read_map_cpp}" \
		" {input.map_encodings_cpp}" \
		" {input.utils_cpp}" \
		" -I {config.cpp_include_dir}" \
		" -I {config.htslib_dir}" \
		" -lhts" \
		" -o {output.bin}"
		
# 2.2 encode genotypes for VCF segments (execute)
rule encode_execute:
	input:
		bin=f"{config.cpp_bin_dir}encode",
		config_file=f"{config.config_file}"
	output:
		encode_log=f"{config.out_dir}encode.log"
	message:
                "Compiling--encoding segments..."
	shell:
		"echo 2. ---ENCODING VCF SEGMENTS---;" \
		"for vcf_f in {config.out_dir}*.vcf; do" \
		"	filename=$(basename $vcf_f);" \
		"	seg_name=${{filename%.*}};" \
		"	echo SEGMENT: $seg_name;" \
		"	./{input.bin} {input.config_file} $vcf_f {config.out_dir}${{seg_name}}.gt {config.out_dir}${{seg_name}}.pos" \
		"	 >> {output.encode_log};" \
		"done;" \
		#"touch output.encode_log;"


# 3.1 build faiss index for encoding segments (compile)
rule faiss_index_compile:
	input:
		encode_log=f"{config.out_dir}encode.log",
		faiss_l2_build_cpp=f"{config.cpp_src_dir}faiss_l2_build.cpp",
		faiss_utils_cpp=f"{config.cpp_src_dir}faiss_utils.cpp"
	output:
		bin=f"{config.cpp_bin_dir}faiss-l2-build"
	message:
		"Compiling--building faiss index for all segments..."
	shell:
		"g++" \
		" {input.faiss_l2_build_cpp}" \
		" {input.faiss_utils_cpp}" \
		" -I {config.cpp_include_dir}" \
		" -I /home/sdp/miniconda3/envs/faiss/include/" \
		" -lfaiss" \
		" -o {output.bin}"	

## 1.1 slice VCF into segments (compile)
#rule slice_VCF_compile:
#	input:
#		interpolate_map=f"{config.data_dir}interpolated.map",
#		main_slice_cpp=f"{config.cpp_src_dir}main_slice.cpp",
#		slice_vcf_cpp=f"{config.cpp_src_dir}slice_vcf.cpp",
#		read_config_cpp=f"{config.cpp_src_dir}read_config.cpp",
#		read_map_cpp=f"{config.cpp_src_dir}read_map.cpp",
#		utils_cpp=f"{config.cpp_src_dir}utils.cpp"
#	output:
#		bin=f"{config.cpp_bin_dir}slice-vcf"
#	message:
#		"Compiling--slice vcf into segments..."
#	shell:
#		"g++" \
#		" {input.main_slice_cpp}" \
#		" {input.slice_vcf_cpp}" \
#		" {input.read_config_cpp}" \
#		" {input.read_map_cpp}" \
#		" {input.utils_cpp}" \
#		" -I {config.cpp_include_dir}" \
#		" -I {config.htslib_dir}" \
#		" -lhts" \
#		" -o {output.bin}"
## 1.2 slice VCF into segments (execute)
#rule slice_VCF_execute:
#	input:
#		bin=f"{config.cpp_bin_dir}slice-vcf",
#		config_file=f"{config.config_file}"
#	output:
#		slice_log=f"{config.out_dir}slice.log"
#	message:
#		"Executing--slice vcf into segments..."
#	shell:
#		"./{input.bin} {input.config_file} > {output.slice_log}"
#
## 2.1 positional encoding vcf segemnts (compile)
#rule pos_encode_vcf_segments_compile:
#	input:
#		slice_log=f"{config.out_dir}slice.log",
#		main_positional_encode_cpp=f"{config.cpp_src_dir}main_positional_encode.cpp",
#		encode_positions_cpp=f"{config.cpp_src_dir}encode_positions.cpp", 
#		map_encodings_cpp=f"{config.cpp_src_dir}map_encodings.cpp",
#		read_config_cpp=f"{config.cpp_src_dir}read_config.cpp",
#		read_map_cpp=f"{config.cpp_src_dir}read_map.cpp",
#		utils_cpp=f"{config.cpp_src_dir}utils.cpp"
#	output:
#		bin=f"{config.cpp_bin_dir}pos-encode"
#	message:
#		"Compiling--positional encode vcf segments..."
#	shell:
#		"g++" \
# 		" {input.main_positional_encode_cpp}" \
#        	" {input.encode_positions_cpp}" \
#        	" {input.map_encodings_cpp}" \
#        	" {input.read_config_cpp}" \
#        	" {input.read_map_cpp}" \
#        	" {input.utils_cpp}" \
#        	" -I {config.cpp_include_dir}" \
#                " -I {config.htslib_dir}" \
#                " -lhts" \
#                " -o {output.bin}"
## 2.2 positional encoding vcf segments (execute)
#rule pos_encode_vcf_segments_execute:	
#	input:
#		bin=f"{config.cpp_bin_dir}pos-encode",
#		config_file=f"{config.config_file}"
#	output:
#		pos_encode_log=f"{config.out_dir}pos-encode.log"
#	message:
#		"Executing--positional encode vcf segments..."
#	shell:
#		"for vcf_f in {config.out_dir}*.vcf; do" \
#                "       filename=$(basename $vcf_f);" \
#                "       seg_name=${{filename%.*}};" \
#		"	echo Positional Encoding: $vcf_f;" \
#		"	num_samples=$(bcftools query -l $vcf_f | wc -l);" \
#                "       ./{input.bin} {input.config_file} $vcf_f {config.out_dir}${{seg_name}}.pos_encoded > {output.pos_encode_log};" \
#                "done;"
#                ##"touch {output.pos_encode_log};"
#		## $bin $config_file $vcf_slice $pos_encode
#
## 3.1 encode vcf segments (compile)
#rule encode_vcf_segments_compile:
#	input:
#		slice_log=f"{config.out_dir}slice.log",
#		main_encode_cpp=f"{config.cpp_src_dir}main_encode.cpp",
#		encode_vcf_cpp=f"{config.cpp_src_dir}encode_vcf.cpp",
#		read_config_cpp=f"{config.cpp_src_dir}read_config.cpp",
#		map_encodings_cpp=f"{config.cpp_src_dir}map_encodings.cpp",
#		utils_cpp=f"{config.cpp_src_dir}utils.cpp"
#	output:
#		bin=f"{config.cpp_bin_dir}encode-vcf"
#	message:
#		"Compiling--encode vcf segments..."
#	shell:
#		"g++" \
#		" {input.main_encode_cpp}" \
#		" {input.encode_vcf_cpp}" \
#		" {input.read_config_cpp}" \
#		" {input.map_encodings_cpp}" \
#		" {input.utils_cpp}" \
#		" -I {config.cpp_include_dir}" \
#		" -I {config.htslib_dir}" \
#		" -lhts" \
#		" -o {output.bin}"
## 3.2 encode vcf segments (execute)
#rule encode_vcf_segments_execute:
#	input:
#		bin=f"{config.cpp_bin_dir}encode-vcf",
#		sample_IDs=f"{config.data_dir}samples_IDs.txt",
#		config_file=f"{config.config_file}"
#	output:
#		encode_log=f"{config.out_dir}encode.log"
#	message:
#		"Executing--encode vcf segments..."
#	shell:
#		"for vcf_f in {config.out_dir}*.vcf; do" \
#		"	filename=$(basename $vcf_f);" \
#		"	seg_name=${{filename%.*}};" \
#		"	./{input.bin} {input.config_file} $vcf_f {config.out_dir}${{seg_name}}.encoded {config.out_dir}${{seg_name}}.position;" \
#		"done;"
#		"touch {output.encode_log};"
#
## 4.1 build FAISS index (compile)
#rule build_faiss_index:
#	input:
#		#slice_log=f"{config.out_dir}slice.log",
#		#pos_encode_log=f"{config.out_dir}pos-encode.log",
#		#encode_log=f"{config.out_dir}encode.log",
#		faiss_l2_build_cpp=f"{config.cpp_src_dir}faiss_l2_build.cpp",
#		faiss_utils_cpp=f"{config.cpp_src_dir}faiss_utils.cpp"
#	output:
#		bin=f"{config.cpp_bin_dir}faiss-l2-build"
#	message:
#		"Compiling--build faiss l2 index..."
#	shell:
#		"g++" \
#		" {input.faiss_l2_build_cpp}" \
#		" {input.faiss_utils_cpp}" \
#		" -I {config.conda_dir}/include" \
#		" -I {config.cpp_include_dir}" \
#		" -I {config.faisslib_dir}" \
#		" -lfaiss" \
#		" -o {output.bin}"
# 
#
## 4 run plink on full vcf
#rule plink:
#	input:
#		vcf=f"{config.vcf_file}"
#	output:
#		plink_done=f"{config.data_dir}plink.log"
#	message:
#		"Running plink --genome on full vcf file"
#	shell:
#		"plink --vcf {input.vcf} --genome --out {config.data_dir}plink"
#
## 5.1 compute euclidean distance for all segments
#rule compute_segment_distance:
#	input:
#		encode_log=f"{config.out_dir}encode.log",
#		query_file=f"{config.query_file}"
#	output:
#		distance_log=f"{config.out_dir}distance.log"
#	message:
#		"Computing Euclidean distance for query against all segments"
#	shell:
#		"for encoded_f in {config.out_dir}*.encoded; do" \
#		"	filename=$(basename $encoded_f);" \
#               "	seg_name=${{filename%.*}};" \
#		"	echo $encoded_f >> {output.distance_log};" \
#		"	python {config.python_dir}distance/compute_segment_distance.py --encoded_file $encoded_f --query_file {input.query_file} > {config.out_dir}${{seg_name}}.dist;" \
#		"done"
## 4.2 aggregate euclidean distance for all segments
#rule aggregate_segment_distance:
#	input:
#		distance_log=f"{config.out_dir}distance.log"
#	output:
#		aggregate_log=f"{config.data_dir}aggregate.log"
#	message:
#		"Aggregating all distance files for all queries"
#	shell:
#		"num_segments=$(ls {config.out_dir}*dist | wc -l)"
#		" && python {config.python_dir}distance/aggregate_segment_distance.py --out_dir {config.out_dir} --ext dist --num_seg $num_segments > {config.data_dir}aggregate.txt"
#		" && touch {config.data_dir}aggregate.log"
#
## 4.2 aggregate euclidean distance for all segments
#rule aggregate_segment_distance:
#	input:
#		distance_log=f"{config.out_dir}distance.log"
#	output:
#		aggregate_log=f"{config.data_dir}aggregate.log"
#	message:
#		"Aggregating all distance files for all queries"
#	shell:
#		"num_segments=$(ls {config.out_dir}*dist | wc -l)"
#		" && python {config.python_dir}distance/aggregate_segment_distance.py --out_dir {config.out_dir} --ext dist --num_seg $num_segments > {config.data_dir}aggregate.txt"
#		" && touch {config.data_dir}aggregate.log"
#
# 5.0 make hap IDs
#rule hap_IDs:
#	input:
#		encode_log=f"{config.out_dir}encode.log"
#	output:
#		hapID_txt=f"{config.data_dir}all_hap_IDs.txt",
#		hapID_log=f"{config.data_dir}hapID.log"
#	message:
#		"generating list of hap IDs..."
#	shell:
#		"for seg_0 in {config.out_dir}*.seg.0.encoded; do" \
#		"	echo $seg_0;" \
#		"	(awk '{{print $1}}' $seg_0) > {config.data_dir}all_hap_IDs.txt;" \
#		"done"
#		" && touch {config.data_dir}hapID.log"
#
## 6.1 faiss (compile)
#rule faiss_compile:
#        input:
#                slice_log=f"{config.out_dir}slice.log",
#                encode_log=f"{config.out_dir}encode.log",
#                faiss_l2_cpp=f"{config.cpp_src_dir}faiss_l2.cpp",
#                build_index_cpp=f"{config.cpp_src_dir}build_index.cpp",
#                read_encodings_cpp=f"{config.cpp_src_dir}read_encodings.cpp",
#                search_index_cpp=f"{config.cpp_src_dir}search_index.cpp",
#                utils_cpp=f"{config.cpp_src_dir}utils.cpp"
#        output:
#                bin=f"{config.cpp_bin_dir}faiss-l2"
#        message:
#                "Compiling--FAISS..."
#        shell:
#                "g++" \
#                " {input.faiss_l2_cpp}" \
#                " {input.build_index_cpp}" \
#                " {input.read_encodings_cpp}" \
#                " {input.search_index_cpp}" \
#                " {input.utils_cpp}" \
#                " -I {config.cpp_include_dir}" \
#                " -I {config.conda_dir}include/" \
#		" -L {config.conda_dir}lib/" \
#                " -lfaiss" \
#                " -o {output.bin}"
#
#
## 6.2 faiss (execute)
#rule faiss_execute:
#	input:
#                bin=f"{config.cpp_bin_dir}faiss-l2",
#		encode_log=f"{config.out_dir}encode.log",
#		database_IDs=f"{config.database_file}",
#		query_IDs=f"{config.query_file}"
#	output:
#		faiss_log=f"{config.out_dir}faiss.log"
#	message:
#		"Executing--FAISS..."
#	shell:
#		"for encoded_f in {config.out_dir}*.encoded; do" \
#		"	filename=$(basename $encoded_f);" \
#		"	seg_name=${{filename%.*}};" \
#               	"	echo $bin ${input.database_IDs} $encoded_f ${input.query_IDs} $encoded_f ${config.k} ${config.delim};" \
#		"done;"
#		"touch {output.faiss_log};"
#                "	seg_name=${{filename%.*}};" \
#		"	echo $encoded_f >> {output.distance_log};" \
#		"	python {config.python_dir}distance/compute_segment_distance.py --encoded_file $encoded_f --query_file {input.query_file} > {config.out_dir}${{seg_name}}.dist;" \
#		"done"
