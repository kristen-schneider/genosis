#ifndef READ_CONFIG_H
#define READ_CONFIG_H

#endif //READ_CONFIG_H

#include <iostream>
#include <fstream>
#include <map>

using namespace std;

map<string,string> get_config_options(string config_file);
map<string,string> get_config_options_old(string config_file);
