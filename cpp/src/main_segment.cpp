#include <cstdlib>
#include <iostream>
#include <fstream>
#include <map>
#include <vector>

#include "read_config.h"
#include "read_map.h"
#include "segment_boundary_map.h"

using namespace std;

int main(int argc, char* argv[]){
    	/*
    	/ takes one config file
    	*/
		
	// read config file
	string configFile = argv[1];   // configuration file will all options
	cout << endl << "1. ---SLICING VCF---" << endl << endl;
	cout << "Loading config options from: " << configFile << "." << endl;
	map<string, string> config_options;
	config_options = get_config_options(configFile);		
	
	// access each option by variable name
	string vcf_file = config_options["vcf_file"];
	string map_file = config_options["interpolated_map"];
	string data_dir = config_options["data_dir"];
	int segment_size = stoi(config_options["segment_size"]);
	
	// report relevant config options
	cout << "CONFIG OPTIONS: " << endl;
	cout << "\t-vcf file: " << vcf_file << endl;
	cout << "\t-map file: " << map_file << endl;
	cout << "\t-data dir: " << data_dir << endl;
	cout << "\t-segment size: " << segment_size << endl << endl;

	// create segment boundary map
	cout << "Making segment boundary map for: " << map_file << " @ " << segment_size << "cm slices." << endl;
	map<int, vector<int>> segment_boundary_map;
	segment_boundary_map = generate_segment_boundary_map(
			map_file,
			segment_size);
	// write segment boundary map
	string segment_boundary_file = data_dir + "segment_boundary.map";
	cout << "Writing segment boundary map to: " << segment_boundary_file << endl;
	write_segment_boundary_map(
			segment_boundary_map,
			segment_boundary_file);


	return 0;
}
