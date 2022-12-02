#ifndef MAP_ENCODINGS_H
#define MAP_ENCODINGS_H

#endif //MAP_ENCODINGS_H

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <map>

using namespace std;

map<string,vector<int>> make_biallelic_encoding_map(string encodingFile);
map<string,int> make_encoding_map(string encodingFile);
