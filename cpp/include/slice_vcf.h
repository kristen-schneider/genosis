#ifndef SLICE_VCF_H
#define SLICE_VCF_H

#endif //SLICE_VCF_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using namespace std;

int slice_main(string map_file, int segment_size, string vcf_file, string out_base_name, string out_dir);
vector<string> read_vcf_header(string vcf_file);
int slice(string vcf_file, vector<string> vcf_header, vector<int> segment_SNP_counts, string base_name, string out_dir);
