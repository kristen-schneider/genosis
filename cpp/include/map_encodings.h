#ifndef MAP_ENCODINGS_H
#define MAP_ENCODINGS_H

#endif //MAP_ENCODINGS_H

#include <iostream>
#include <fstream>
#include <map>

using namespace std;

map<string,int> make_encoding_map(string encodingFile);
map<string,vector<int>> make_biallelic_encoding_map(string encodingFile);
void split_line(const string &s, char delim, vector<string> &elems);
