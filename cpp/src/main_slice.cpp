#include <iostream>
#include <fstream>
#include <map>
#include <vector>

#include "map_encodings.h"
#include "read_config.h"
#include "slice_vcf.h"

using namespace std;

int main(int argc, char* argv[]){
    	/*
    	/ takes one config file
    	*/
		
	// read config file
	cout << "Loading Config Options..." << endl;
	string configFile = argv[1];   // configuration file will all options
	map<string, string> config_options;
	config_options = get_config_options(configFile);		
	// access each option by variable name
	string vcf_file = config_options["vcf_file"];
	string map_file = config_options["map_file"];
	string out_dir = config_options["out_dir"];
    	string out_base_name = config_options["out_base_name"];
	int segment_size = stoi(config_options["segment_size"]);
	
	// slice vcf into segments
	cout << "Slicing VCF..." << endl;
	int num_segments = slice_main(map_file, segment_size, vcf_file, out_base_name, out_dir);
	cout << "Wrote " << num_segments + 1 << " slices." << endl; // zero-index
	
	return 0;
}
