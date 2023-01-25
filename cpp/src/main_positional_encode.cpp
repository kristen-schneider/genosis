#include <cstdlib>
#include <iostream>
#include <fstream>
#include <map>
#include <vector>

#include "encode_positions.h"
#include "read_map.h"
#include "read_config.h"

using namespace std;

int main(int argc, char* argv[]){
    	/*
    	/ takes one config file
    	*/
		
	// commandline arguments
	string configFile = argv[1];   		// configuration file will all options
	string vcf_slice_file = argv[2];	// input vcf (slice) file
	//string output_encoding_file = argv[3];	// output encoding (slice) file
	string output_position_file = argv[3];	// output positional encoding (slice) file
	
	// options from config file
	cout << "Loading Config Options..." << endl;
	map<string, string> config_options;
	config_options = get_config_options(configFile);		
	// access each option by variable name
	string map_file = config_options["map_file"];
	string encoding_file = config_options["encoding_file"];
	string sample_IDs_file = config_options["sample_IDs_file"];
	string out_dir = config_options["out_dir"];
    	string out_base_name = config_options["out_base_name"];
	
	// make positional  map
	encode_positions(map_file);

	// encode single vcf
	//encode_vcf(sample_IDs_file, vcf_slice_file, encoding_map, output_encoding_file, output_position_file);	
	
	return 0;
}
