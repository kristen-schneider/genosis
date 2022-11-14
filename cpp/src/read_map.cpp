#include <iostream>
#include <fstream>
#include <string>

#include "read_map.h"
#include "utils.h"

using namespace std;

/*
 * Open a map file and read each SNP record
 * count SNPs until XcM is reached
 * return a vector of SNP lengths for each XcM slice
*/
vector<int> read_map_file(string map_file, float slice_size){
	float max_cm = 1.0;
	int snp_count = 0;
	int slice_count = 0;

	vector<int> slice_SNP_counts;
	
	ifstream map_file_stream;
        map_file_stream.open(map_file);
        if (!map_file_stream.is_open()){
                cout << "FAILED TO OPEN: " << map_file << endl;
                exit(1);
        }
	cout << "...Reading map file..." << endl;
	string line;
	while (getline (map_file_stream, line)){
		snp_count ++;
		int cm_index = 2;
		vector<string> single_SNP;
        	split_line(line, ' ', single_SNP);
		float record_cm = stof(single_SNP[cm_index]);
		if (record_cm >= max_cm){
			slice_SNP_counts.push_back(snp_count);
			snp_count = 0;
			max_cm = record_cm + slice_size;
			slice_count ++;
		}
	}
        slice_SNP_counts.push_back(snp_count);
        snp_count = 0;
	cout << "...Done reading map file." << endl;
	return slice_SNP_counts;
}

/*
 * Splits a string by some delimiter
 
void split_line(const string &s, char delim, vector<string> &elems){
    stringstream ss;
    ss.str(s);
    string item;
    while (getline(ss, item, delim)) {
	    elems.push_back(item);
    }
}*/
