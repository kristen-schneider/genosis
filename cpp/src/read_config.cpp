#include <iostream>
#include <fstream>
#include <map>

#include "read_config.h"
#include "utils.h"

using namespace std;

/*
 * Reads configuration file and stores options in a map
 *
 * @param config_file: (string) path to config file
 * @returns config_options: (map<string, string>) 
 * 		key: option keyword
 * 		value: user-specified option
 */
map<string, string> get_config_options(string config_file){
	// to return
	map <string, string> config_options;

	// open config file and check success
	ifstream config_file_stream;
	config_file_stream.open(config_file);
	if (!config_file_stream.is_open()){
		cout << "FAILED TO OPEN: " << config_file << endl;
		exit(EXIT_FAILURE);
	}
		
	char delim = ':';
	string line;
	string key;
	string user_specified_option;

	while (getline(config_file_stream, line)){
		if (line.find(delim) !=  std::string::npos){
			vector<string> _line;
			split_line(line, delim, _line);
			key = _line[0];
			continue;
		}else{
			user_specified_option = line;
		}
		pair<string, string> key_option_pair (key, user_specified_option);
		config_options.insert(key_option_pair);
	}
	return config_options;
}


map<string, string> get_config_options_old(string config_file){

	// a map which stores configuration options
	map <string, string> config_options;
	
	// delimiter for config file
	string delim = ":";

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
