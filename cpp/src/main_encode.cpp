#include <cstdlib>
#include <iostream>
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

        	
	string full_vcf_segment_file = vcf_segment_file;    // store full path to vcf segment file
                                                            // the original input string is manipulated below.
	
        // 1. get output file name
        // split the file path and get the file name
	vector<string> file_path;
	string path_delim = "/";
	size_t pos = 0;
	string path_token;
	while ((pos = vcf_segment_file.find(path_delim)) != std::string::npos){
		path_token = vcf_segment_file.substr(0, pos);
		file_path.push_back(path_token);
		vcf_segment_file.erase(0, pos + path_delim.length());
	}
	path_token = vcf_segment_file.substr(0, pos);
	file_path.push_back(path_token);
	// split the file name and get the chrm idx and seg idx
	string file_name = file_path.back();
	vector<string> file_name_vec;
	string file_delim = ".";
	pos = 0;
	string file_token;
	while ((pos = file_name.find(file_delim)) != std::string::npos){
                file_token = file_name.substr(0, pos);
                file_name_vec.push_back(file_token);
                file_name.erase(0, pos + file_delim.length());
        }
        file_token = file_name.substr(0, pos);
        file_name_vec.push_back(file_token);
	int chrm_idx = stoi(file_name_vec[0].erase(0, 4));      // chromosome idx
	int segment_idx = stoi(file_name_vec[1].erase(0, 7));   // segment idx
	// name output files
	string output_gt_file = encoding_dir + "chrm" + to_string(chrm_idx) + ".segment" + to_string(segment_idx) + ".gt";
	string output_pos_file = encoding_dir + "chrm" + to_string(chrm_idx) + ".segment" + to_string(segment_idx) + ".pos";

	
	// 2. make gt encoding map
	// format: (gt (string): encoding <int, int>)
	map<string, vector<int>> encoding_def_map = map_gt_encoding(encoding_def_file);	
	
        // 3. make cm positional encoding map
	// format: (chrm (int): bp_cm <int, float>)
	map<int, map<int, float>> chrm_bp_cm_map = map_bp_cm(interpolated_map);
	
	// 4. encode single vcf
	encode_vectors(chrm_idx,
			full_vcf_segment_file,
			sample_IDs_file,
			encoding_def_map,
			chrm_bp_cm_map,
			output_gt_file,
			output_pos_file);	
	
	return 0;
}
