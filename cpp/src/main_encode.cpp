#include <cstdlib>
#include <iostream>
#include <fstream>
#include <map>
#include <vector>

#include "encode_gt.h"
#include "map_encodings.h"
#include "read_config.h"

using namespace std;

int main(int argc, char* argv[]){
    	/*
    	/ takes one config file
    	*/
		
	// read config file
	string configFile = argv[1];   		// configuration file will all options
	string vcf_slice_file = argv[2];	// input vcf (slice) file
	string output_encoding_file = argv[3];	// output encoding (slice) file
	string output_position_file = argv[4];	// output positional encoding (slice) file
	map<string, string> config_options;
	config_options = get_config_options(configFile);		
	
	// access each option by variable name
	string map_file = config_options["interpolated_map"];
	string encoding_file = config_options["encoding_file"];
	string sample_IDs_file = config_options["sample_IDs_file"];
	string out_dir = config_options["out_dir"];
    	string out_base_name = config_options["out_base_name"];

	// make encoding map
	cout << "...Loading encoding map." << endl;
	map<string, vector<int>> encoding_map = make_biallelic_encoding_map(encoding_file);
	cout << "...Done loading encoding map." << endl;

	// encode single vcf
	encode_gt_vectors(sample_IDs_file,
			vcf_slice_file,
			encoding_map,
			output_encoding_file);	
	
	return 0;
}
