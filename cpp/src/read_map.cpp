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
			bp_start_end.clear();
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
 * read map file, map bp pos to cm pos
 */
map<int, float> make_bp_cm_map(string map_file){
	map<int, float> bp_cm_map;
	ifstream map_file_stream(map_file);
	string line;
        while (getline (map_file_stream, line)){

                // column with cm and bp data
                float cm_col = 2;
                int bp_col = 3;

                // read and split line
                vector<string> map_line;
                split_line(line, ' ', map_line);
		// extract bp and cm values
		int bp_value = stoi(map_line[bp_col]);
		float cm_value = stof(map_line[cm_col]);
		// add pair to map
        	bp_cm_map[bp_value] = cm_value;
	}


    return bp_cm_map;
}
