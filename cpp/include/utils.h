#ifndef UTILS_H
#define UTILS_H

#endif //UTILS_H

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

vector<vector<int>> transpose(vector<vector<int>> &vmf);
vector<string> get_sample_IDs(string sample_IDs_file);
int get_num_samples(bcf_hdr_t *vcf_header);
void get_dimensions(string in_file, int* dimensions, char delim);
