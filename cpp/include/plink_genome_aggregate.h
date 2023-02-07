#ifndef PLINK_GENOME_AGGREGATE_H
#define PLINK_GENOME_AGGREGATE_H

#endif //PLINK_GENOME_AGGREGATE_H

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <map>
#include <vector>


using namespace std;

map<string, map<string, float>> main_plink_aggregate(string plink_file, char delim);
vector<pair<string, float>> sort_map(map<string, float> str_flt_map);
map<string, vector<pair<string, float>>> sort_full_map(map<string, map<string, float>> str_str_flt_map);
map<string, vector<pair<string, float>>> return_top_k(map<string, vector<pair<string, float>>> str_vec_map, vector<string> queries, int k);
void write_top_k(map<string, vector<pair<string, float>>> top_k_map, int k, string output_file);
void split_line(const string &s, char delim, vector<string> &elems);
vector<string> read_query_file(string query_file);
