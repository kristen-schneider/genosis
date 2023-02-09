#include <cstdlib>
#include <iostream>
#include <fstream>
#include <map>
#include <vector>

#include "read_config.h"
#include "read_map.h"
#include "slice_vcf.h"

using namespace std;

int main(int argc, char* argv[]){
    	/*
    	/ takes one config file
    	*/
		
	// read config file
	string configFile = argv[1];   // configuration file will all options
	cout << endl << "1. ---SLICING VCF---" << endl << endl;
	cout << "Loading config options from: " << configFile << "..." << endl;
	map<string, string> config_options;
	config_options = get_config_options(configFile);		
	
	// access each option by variable name
	string vcf_file = config_options["vcf_file"];
	string map_file = config_options["interpolated_map"];
	string out_dir = config_options["out_dir"];
    	string out_base_name = config_options["out_base_name"];
	int segment_size = stoi(config_options["segment_size"]);
	
	// report relevant config options
	cout << "CONFIG OPTIONS: " << endl;
	cout << "\t-vcf file: " << vcf_file << endl;
	cout << "\t-map file: " << map_file << endl;
	cout << "\t-out dir: " << out_dir << endl;
	cout << "\t-segment size: " << segment_size << endl << endl;

	// slice vcf into segments
	cout << "Slicing VCF: " << vcf_file << endl;
	
	int num_segments = slice_main(map_file,
			segment_size,
			vcf_file,
			out_base_name,
			out_dir);
	cout << "Wrote " << num_segments + 1 << " slices to: " << out_dir << endl; // zero-index
	
	return 0;
}
