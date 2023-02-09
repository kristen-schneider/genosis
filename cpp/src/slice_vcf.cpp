#include <iostream>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>
#include <cstdlib>

#include "read_map.h"
#include "slice_vcf.h"
#include "utils.h"

using namespace std;

/*
 * Slice vcf file into segments of specified cm length
 */	
int slice_main(string map_file,
		int segment_size,
		string vcf_file,
		string out_base_name,
		string out_dir){

	// slice full chromosome VCF file into smaller slices
	cout << "...preparing to slice VCF into segments of size: " << segment_size << "cM." << endl;
 
	// map segment index to bp start, bp end
	map<int, vector<int>> segment_boundary_map;
	segment_boundary_map = make_segment_boundary_map(map_file,
			segment_size);
	int num_segments = slice(vcf_file,
			segment_boundary_map,
			out_base_name,
			out_dir);
		
	cout << "...Done reading VCF file." << endl;
	cout << "Done slicing." << endl;
	return num_segments;
	
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
