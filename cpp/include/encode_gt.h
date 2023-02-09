#ifndef ENCODEVCF_H
#define ENCODEVCF_H

#endif //ENCODEVCF_H

#include <cstdlib>
#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <htslib/vcfutils.h>
#include <htslib/kstring.h>
#include <htslib/kseq.h>
#include <htslib/synced_bcf_reader.h>
#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <vector>

using namespace std;

void encode_gt_vectors(string sample_IDs_file,
	       	string input_vcf_file, 
		map<string, vector<int>> encoding_map, 
		string output_encoding_file);

void write_SMF(vector<string> all_sample_IDs, 
		vector<vector<int>> smf, 
		string output_encoding_file);

//void write_SMF_haplotype(vector<string> all_sample_IDs, vector<vector<int>> smf, string output_encoding_file);
void write_positional_encoding(vector<int> all_positions,
                vector<string> all_sample_IDs,
                vector<vector<int>> smf, 
		string output_position_file);
