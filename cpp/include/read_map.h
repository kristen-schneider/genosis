#ifndef READ_MAP_H
#define READ_MAP_H

#endif //READ_MAP_H

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using namespace std;

map<int,vector<int>> make_segment_boundary_map(string map_file, int slice_size);
vector<int> read_map_file(string map_file, float slice_size);
map<int, float> make_bp_cm_map(string map_file);
