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

int encode_main();
void encode_vcf(string sample_IDs_file, string input_vcf_file, map<string, vector<int>> encoding_map, string output_encoding_file);
void write_SMF(vector<string> all_sample_IDs, vector<vector<int>> smf, string output_encoding_file);
void write_SMF_haplotype(vector<string> all_sample_IDs, vector<vector<int>> smf, string output_encoding_file);
int get_num_samples(bcf_hdr_t *vcf_header);
const char **get_sequence_names(bcf_hdr_t *vcf_header);
