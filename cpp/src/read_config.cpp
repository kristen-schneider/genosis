#include <iostream>
#include <fstream>
#include <map>

#include "read_config.h"

using namespace std;

map<string, string> get_config_options(string config_file){
	/*
	 * Opens a configuration file and stores options in map.
	 */

	// a map which stores configuration options
	map <string, string> config_options;
	
	// delimiter for config file
	string delim = "=";

	// open file and check success
	ifstream config_file_stream;
	config_file_stream.open(config_file);
	if (!config_file_stream.is_open()){
		cout << "FAILED TO OPEN: " << config_file << endl;
	}

	// read file
	else{
		string line;
		string key;
		string option;

		// split by delim
		auto start = 0U;
		auto end = 0U;
		while(getline(config_file_stream, line)){
			// get key and option from file
			end = line.find(delim, start);
			key = line.substr(start, end);	// key
			start = end + delim.length();	
			option = line.substr(start);	// option
			start = 0U;

			// add key and option to map
            		pair<string, string> p (key, option);
            		config_options.insert(p);
		}
	}
		
	return config_options;
}
