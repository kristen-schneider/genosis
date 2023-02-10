#include <cstdlib>
#include <iostream>
#include <fstream>
#include <map>
#include <vector>

#include "encode_segment.h"
#include "map_encodings.h"
#include "read_config.h"
#include "read_map.h"

using namespace std;

int main(int argc, char* argv[]){
    	/*
    	/ takes one config file
    	*/
		
	// read config file
	string configFile = argv[1];   		// configuration file will all options
	string vcf_slice_file = argv[2];	// input vcf (slice) file
	string output_gt_file = argv[3];	// output encoding (slice) file
	string output_pos_file = argv[4];	// output positional encoding (slice) file
	map<string, string> config_options;
	config_options = get_config_options(configFile);		
	
	// access each option by variable name
	string map_file = config_options["interpolated_map"];
	string encoding_file = config_options["encoding_file"];
	string sample_IDs_file = config_options["sample_IDs_file"];
	string out_dir = config_options["out_dir"];
    	string out_base_name = config_options["out_base_name"];

	// make nominal gt encoding map
	cout << "...Loading gt encoding map." << endl;
	map<string, vector<int>> encoding_map = make_biallelic_encoding_map(encoding_file);
	cout << "...Done loading gt encoding map." << endl;
	// make cm positional encoding map
	cout << "...Loading pos encoding map." << endl;
	map<int, float> bp_cm_map;
	bp_cm_map = make_bp_cm_map(map_file);
	cout << "...Done loading pos encoding map." << endl;
	
	// encode single vcf
	encode_vectors(sample_IDs_file,
			vcf_slice_file,
			encoding_map,
			bp_cm_map,
			output_gt_file,
			output_pos_file);	
	
	
	return 0;
}
