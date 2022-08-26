#ifndef SLICE_VCF_H
#define SLICE_VCF_H

#endif //SLICE_VCF_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using namespace std;

vector<string> read_vcf_header(string vcf_file);
int slice(string vcf_file, vector<string> vcf_header, int segment_size, string base_name, string out_dir);;

