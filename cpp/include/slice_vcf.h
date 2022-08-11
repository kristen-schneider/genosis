#ifndef SLICE_VCF_H
#define SLICE_VCF_H

#endif //SLICE_VCF_H

#include <iostream>
#include <fstream>
#include <string>

using namespace std;

void write_vcf_header(string vcf_file, string out_file);
int slice(string vcf_file, int segment_size, string base_name, string out_dir);;
