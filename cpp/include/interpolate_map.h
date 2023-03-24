#ifndef INTERPOLATE_MAP_H
#define INTERPOLATE_MAP_H

#endif //INTERPOLATE_MAP_H

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iterator>
#include <map>

using namespace std;

// map of chromosomes to vectors of basebpairs
map<string, vector<int>> make_chr_bp_map(
        string vcf_bps);
// map (basepair: centimorgan) map for each chromosome
map<int, vector<tuple<int, float>>> make_chr_cm_map(
        string map_file);
// interpolate map and write new output file
void interpolate_map(
        map<string, vector<int>> chr_bp_map,
        map<int, vector<tuple<int, float>>> chrm_cm_map,
        string outfile);
