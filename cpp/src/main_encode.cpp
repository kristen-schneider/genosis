#include <cstdlib>
#include <iostream>
#include <sstream>
#include <fstream>
#include <map>
#include <vector>

#include "encode_segment.h"
#include "map_encodings.h"
#include "read_map.h"
#include "utils.h"

using namespace std;

int main(int argc, char* argv[]){
	
        // 0. get input options
	string vcf_segment_file = argv[1];	// input vcf (segment) file to encode
	string sample_IDs_file = argv[2];	// sample IDs to pair with encodings
	string encoding_def_file = argv[3];	// type of encodings to use (encoding def)
	string interpolated_map= argv[4];	// interpolated map file
	string encoding_dir = argv[5];		// directory to which encoding files should be written out
	   
	// 1. get output file name
	stringstream full_vcf_segment_path(vcf_segment_file);	// store full path to vcf segment file
	string vcf_segment_basename;
	while(getline(full_vcf_segment_path, vcf_segment_basename, '/')){}	// get the last token in the path
	// get chromosome and segment from vcf_segment_basename
	stringstream vcf_segment_basename_ss(vcf_segment_basename);
	// format: chrm1.segment0.vcf.gz
	vector<string> vcf_segment_basename_vec;
	string vcf_segment_token;
	while(getline(vcf_segment_basename_ss, vcf_segment_token, '.')){
		vcf_segment_basename_vec.push_back(vcf_segment_token);
	}

	int chrm_idx = stoi(vcf_segment_basename_vec[0].erase(0, 4));	// chromosome idx
	int segment_idx = stoi(vcf_segment_basename_vec[1].erase(0, 7));	// segment idx

	// name output files
	string output_gt_file = encoding_dir + "chrm" + to_string(chrm_idx) + ".segment" + to_string(segment_idx) + ".gt";
	string output_pos_file = encoding_dir + "chrm" + to_string(chrm_idx) + ".segment" + to_string(segment_idx) + ".pos";
	
	// 2. make gt encoding map
	// format: (gt (string): encoding <int, int>)
	map<string, vector<int> > encoding_def_map = map_gt_encoding(encoding_def_file);	
	
    // 3. make cm positional encoding map
	// format: (chrm (int): bp_cm <int, float>)
	map<int, map<int, float> > chrm_bp_cm_map = map_bp_cm(interpolated_map);
	
	// 4. encode single vcf
	encode_vectors(chrm_idx,
			vcf_segment_file,
			sample_IDs_file,
			encoding_def_map,
			chrm_bp_cm_map,
			output_gt_file,
			output_pos_file);	
	
	return 0;
}
