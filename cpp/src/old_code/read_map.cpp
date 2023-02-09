#include <cstdlib>
#include <iostream>
#include <fstream>
#include <map>
#include <string>

#include "read_map.h"
#include "utils.h"

using namespace std;

/*
 * Open a map file and real each SNP record
 * return a map of bp start and end bp positions
 * cm_idx: <start_bp, end_bp>
*/
map<int,vector<int>> make_segment_boundary_map(string map_file,
		int slice_size){

	map<int, vector<int>> segment_boundary_map; 	// cm_index: <start_bp, end_bp>
	vector<int> bp_start_end; 			// <start_bp, end_bp>
	
	int seg_index = 0;
	int bp_start = 0;
	int bp_end = -1;
	
	float curr_cm;
	int curr_bp;
	float curr_slice_size;

	// open map file
	ifstream map_file_stream;
        map_file_stream.open(map_file);
        if (!map_file_stream.is_open()){
                cout << "FAILED TO OPEN MAP FILE: " << map_file << endl;
                exit(1);
        }
	
	// read map file
        cout << "...Reading map file." << endl;
        cout << "......Creating segment boundary map (segment index: start_bp, end_bp)" << endl;
        string line;
        while (getline (map_file_stream, line)){
                
		// column with cm and bp data
		float cm_col = 2;
		int bp_col = 3;

		// read and split line
               	vector<string> map_line;
                split_line(line, ' ', map_line);
                
		// read values from line
		curr_cm = stof(map_line[cm_col]);		// current cm data
		curr_bp = stoi(map_line[bp_col]);		// current bp data
		curr_slice_size = curr_cm - seg_index;		// current slice size size
                
		// when a full slice is found...
		if (curr_slice_size >= slice_size){
			// populate start and end vector
			bp_start_end.push_back(bp_start);
			bp_end = curr_bp;
			bp_start = bp_end;
			bp_start_end.push_back(bp_end);
		
			// fill out the map 
			segment_boundary_map[seg_index] = bp_start_end;
			seg_index += 1;
                }
        }
	// last segment
	// populate start and end vector
        bp_start_end.push_back(bp_start);
        bp_end = curr_bp;
        bp_start = bp_end;
        bp_start_end.push_back(bp_end);

        // fill out the map 
        segment_boundary_map[seg_index] = bp_start_end;
        seg_index += 1;
	cout << "......counted " << segment_boundary_map.size() << " slices." << endl;
	cout << "...Done reading map file." << endl;
        return segment_boundary_map;
}


/*
 * reads a map file and returns a < pos, cm > map
 * which maps a bp pos to a cm pos from the map file
 */
map<int, float> make_bp_segment_boundary_map(string map_file){

        // map data structure <bp, cm>
        map<int, float> bp_segment_boundary_map;

	// open map file
	ifstream map_file_stream;
        map_file_stream.open(map_file);
        if (!map_file_stream.is_open()){
                cout << "FAILED TO OPEN MAP FILE: " << map_file << endl;
                exit(1);
        }
	
	// read map file
        cout << "...Reading map file..." << endl;
        cout << "...mapping positional encodings..." << endl;
        string line;
        while (getline (map_file_stream, line)){
                // column with cm and bp data
		float cm_col = 2;
		int bp_col = 3;

		// read and split line
               	vector<string> map_line;
                split_line(line, ' ', map_line);
		float curr_cm = stof(map_line[cm_col]);	// current cm data
		int curr_bp = stoi(map_line[bp_col]);	// current bp data

		bp_segment_boundary_map[curr_bp] = curr_cm;
	}
	return bp_segment_boundary_map;
}


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
                cout << "FAILED TO OPEN MAP FILE: " << map_file << endl;
                exit(1);
        }
	cout << "......Reading map file..." << endl;
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
	return slice_SNP_counts;
}
