#ifndef ENCODEVCF_H
#define ENCODEVCF_H

#endif //ENCODEVCF_H

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

void write_all_segments(int num_segments, map<string, int> encoding_map, string output_dir, string output_base_name);
void write_encoded_vcf(string input_vcf_file, map<string, int> encoding_map, string output_encoding_file);
int get_num_samples(bcf_hdr_t *vcf_header);
const char **get_sequence_names(bcf_hdr_t *vcf_header);
vector<vector<int>> transpose(vector<vector<int>> &vmf);
void write_SMF(vector<vector<int>> smf, string output_encoding_file);
