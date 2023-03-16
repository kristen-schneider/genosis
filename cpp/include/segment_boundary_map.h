#ifndef SEGMENT_BOUNDARY_MAP_H
#define SEGMENT_BOUNDARY_MAP_H

#endif //SEGMENT_BOUNDARY_MAP_H

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iterator>
#include <map>
#include <tuple>

using namespace std;

map<int, vector<tuple<int, int>>> make_slice_boundary_map(
        string interpolated_map,
        float slice_size);
void write_slice_boundary_map(
        map<int, vector<tuple<int, int>>> slice_boundary_map,
        string outfile);
