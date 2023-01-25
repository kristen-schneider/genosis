#ifndef UTILS_H
#define UTILS_H

#endif //UTILS_H

#include <cstdlib>
#include <iostream>
#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <htslib/vcfutils.h>
#include <htslib/kstring.h>
#include <htslib/kseq.h>
#include <htslib/synced_bcf_reader.h>
#include <string>
#include <vector>
#include <map>


using namespace std;

vector<vector<int>> transpose_int(vector<vector<int>> &vmf);
vector<vector<float>> transpose_float(vector<vector<float>> &vmf);
vector<string> get_sample_IDs(string sample_IDs_file);
int get_num_samples(bcf_hdr_t *vcf_header);
int count_num_samples(string in_file);
void split_line(const string &s, char delim, vector<string> &elems);
int count_length_input_vector(string in_file, char delim);
const char **get_sequence_names(bcf_hdr_t *vcf_header);
