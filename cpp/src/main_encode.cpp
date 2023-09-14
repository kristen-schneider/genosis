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
	// get input options
	string vcf_segment_file = argv[1];	// input vcf (segment) file to encode
	string sample_IDs_file = argv[2];	// sample IDs to pair with encodings
	string encoding_file = argv[3];		// encodings to use
	string interpolated_map= argv[4];	// where interpolated map file should exist
	string encoding_dir = argv[5];		// where the encoding files should be written out
	//string output_gt_file = argv[6];	// output gt encoded (segment) file
	//string output_pos_file = argv[7];	// output pos encoded (segment) file
	//string output_af_file = argv[8];	// output allele frequency encoding (segment) file
	
	string full_vcf_segment_file = vcf_segment_file;	
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
	int chrm_idx = stoi(file_name_vec[0].erase(0, 4));
	int segment_idx = stoi(file_name_vec[1].erase(0, 7));

	// name output files
	string output_gt_file = encoding_dir + "chrm" + to_string(chrm_idx) + ".segment" + to_string(segment_idx) + ".gt";
	string output_pos_file = encoding_dir + "chrm" + to_string(chrm_idx) + ".segment" + to_string(segment_idx) + ".pos";

	
	// make gt encoding map (gt (string): encoding <int, int>)
	//cout << "...Loading gt encoding map." << endl;
	map<string, vector<int>> encoding_map = map_gt_encoding(encoding_file);
	//cout << "...Done loading gt encoding map." << endl;
	
	// make cm positional encoding map
	//cout << "...Loading pos encoding map." << endl;
	map<int, map<int, float>> bp_cm_map = map_bp_cm(interpolated_map);
	//cout << "...Done loading pos encoding map." << endl;
	
	cout << vcf_segment_file << " " << full_vcf_segment_file << endl;
	
	
	// encode single vcf
	encode_vectors(chrm_idx,
			full_vcf_segment_file,
			sample_IDs_file,
			encoding_map,
			bp_cm_map,
			output_gt_file,
			output_pos_file);	
	
	return 0;
}
