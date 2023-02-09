#include <iostream>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>
#include <cstdlib>

#include "segment_boundary_map.h"
#include "read_map.h"
#include "utils.h"

using namespace std;

/*
 * creating segment boundary map with format
 * segment index: <bp_start, bp_end>
 */	
map<int, vector<int>> make_segment_boundary_map(string map_file,
		int segment_size,
		string out_dir){
	
	// map segment index to bp start, bp end
	map<int, vector<int>> segment_boundary_map;
	segment_boundary_map = make_segment_boundary_map(map_file,
			segment_size);
	/*
	int num_segments = slice(vcf_file,
			segment_boundary_map,
			out_base_name,
			out_dir);
		
	return num_segments;
	*/
	return segment_boundary_map;
}

/*
 * write segment boundary map to file 
 * segment index, bp_start, bp_end 
 */
void write_segment_boundary_map(map<int, vector<int>> segment_boundary_map, 
		string segment_boundary_file){
	// open out file
	ofstream segment_boundary_file_stream;
	
	segment_boundary_file_stream.open(segment_boundary_file);
	
	for (auto it = segment_boundary_map.begin(); it != segment_boundary_map.end(); ++it){
		segment_boundary_file_stream << it->first << " ";
		for (int i = 0; i < it->second.size(); i++){
			segment_boundary_file_stream << it->second[i] << " ";
		}
		segment_boundary_file_stream << endl;
	}
}


/*
 * Open full VCF file, ignore header, 
 * write a one slice to slice file
 * return number of slices
 */
int slice(string vcf_file, 
		map<int, vector<int>> cm_map, 
		string base_name,
		string out_dir){

	int num_segments = -1;

	// get VCF header from file
	const char *vcf_file_cc = vcf_file.c_str();
	htsFile *vcf_stream = bcf_open(vcf_file_cc, "r");
	bcf_hdr_t *vcf_header = bcf_hdr_read(vcf_stream);
	if (vcf_header == NULL) {
		// if VCF header is empty, print to screen
		cout << "...read empty  VCF header." << endl;
		//throw runtime_error("Unable to read header.");
		//exit(1);
	}

	return num_segments;
}
