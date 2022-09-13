#ifndef SLICE_VCF_H
#define SLICE_VCF_H

#endif //SLICE_VCF_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using namespace std;

vector<string> read_vcf_header(string vcf_file);
int slice(string vcf_file, vector<string> vcf_header, vector<int> segment_SNP_counts, string base_name, string out_dir);
vector<int> read_map_file(string map_file, float slice_size);
void split_line(const string &s, char delim, vector<string> &elems);
