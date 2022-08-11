#include <iostream>
#include <fstream>
#include <map>
#include <vector>

#include "encode_vcf.h"
#include "map_encodings.h"
#include "read_config.h"
#include "slice_vcf.h"

using namespace std;

int main(int argc, char* argv[]){
    	/*
    	/ takes one config file
    	*/
		
	// read config file
	cout << "Reading Config File..." << endl;
	string configFile = argv[1];   // configuration file will all options
	map<string, string> config_options;
	config_options = get_config_options(configFile);
	
	// access each option by variable name
	string encoding_file = config_options["encoding_file"];
	string vcf_file = config_options["vcf_file"];
	int segment_size = stoi(config_options["segment_size"]);
	string out_dir = config_options["out_dir"];
    	string output_base_name = config_options["out_base_name"];
	
	// make encoding map
	map<string, int> encoding_map = make_encoding_map(encoding_file);

	// slice vcf into segments
	cout << "Slicing VCF..." << endl;
	int num_segments = slice(vcf_file, segment_size, output_base_name, out_dir);


	// encoded vcf segments 
	cout << "Encoding VCF..." << endl;
	write_all_segments(num_segments, encoding_map, out_dir, output_base_name);

	return 0;
}
