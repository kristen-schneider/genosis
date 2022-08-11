#configfile: "config.yaml"

rule all:
	input:
		"data/segments/seg_5000/segments.encoding.done"

rule segment_vcf_COMPILE:
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
		

rule segment_vcf_RUN:
	input:
		bin="cpp/bin/segment",
		config_file="cpp/configs/segment_config"
	output:
		done="data/segments/seg_5000/segments.encoding.done"
	message: 
		"Executing--slice vcf into segments and encode"
	shell:
		"./{input.bin} {input.config_file} && touch {output.done}"
