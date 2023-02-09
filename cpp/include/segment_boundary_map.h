#ifndef SLICE_VCF_H
#define SLICE_VCF_H

#endif //SLICE_VCF_H

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <vector>

using namespace std;

map<int, vector<int>> make_segment_boundary_map(string map_file,
                int segment_size,
                string out_dir);
void write_segment_boundary_map(map<int, vector<int>> segment_boundary_map,
                string segment_boundary_file);

//int slice_main(string map_file, int segment_size, string vcf_file, string out_base_name, string out_dir);
//int slice(string vcf_file, map<int, vector<int>> cm_map, string base_name, string out_dir);
