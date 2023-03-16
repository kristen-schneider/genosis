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
	int chrm_idx = stoi(argv[1]);
	string vcf_segment_file = argv[2];	// input vcf (segment) file to encode
	string sample_IDs_file = argv[3];	// sample IDs to pair with encodings
	string encoding_file = argv[4];		// encodings to use
	string interpolated_map= argv[5];	// where interpolated map file should exist
	string output_gt_file = argv[6];	// output gt encoded (segment) file
	string output_pos_file = argv[7];	// output pos encoded (segment) file
	string output_af_file = argv[8];	// output allele frequency encoding (segment) file
	
	// make gt encoding map (gt (string): encoding <int, int>)
	cout << "...Loading gt encoding map." << endl;
	map<string, vector<int>> encoding_map = map_gt_encoding(encoding_file);
	cout << "...Done loading gt encoding map." << endl;
	
	// make cm positional encoding map
	cout << "...Loading pos encoding map." << endl;
	map<int, map<int, float>> bp_cm_map = map_bp_cm(interpolated_map);
	cout << "...Done loading pos encoding map." << endl;
	
	// make af encoding map
	
	
	// encode single vcf
	encode_vectors(chrm_idx,
			vcf_segment_file,
			sample_IDs_file,
			encoding_map,
			bp_cm_map,
			output_gt_file,
			output_pos_file,
			output_af_file);	
	
	
	return 0;
}
