#include <cstdlib>
#include <faiss/IndexFlat.h>
#include <faiss/IndexHNSW.h>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>

#include "build_index.h"


/**
 * Copyright (c) Facebook, Inc. and its affiliates.
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


// 64-bit int
using idx_t = faiss::Index::idx_t;
using namespace std;

faiss::IndexFlatL2 build_l2_index(string database_IDs, string database_data, char delim, int num_elements){
	// make map of sample ID to sample encoding (embedding)
	map<string, float*> ID_database_map = make_ID_data_map(database_data, delim, num_elements);
	
	faiss::IndexFlatL2 index(num_elements);
	//int i = 0;
	//float* sample_arr = new float[num_elements];

	// read through databse file (list of samples in database)
	// and pull those vectors from the map to add to index
	ifstream db_file_stream;
	db_file_stream.open(database_IDs);
	if (!db_file_stream.is_open()){
		cout << "Failed to open: " << database_IDs << endl;
		exit(1);
	}
	string line;
	while (getline(db_file_stream, line)){
		// pull sample vector from map
		float* sample_vector = ID_database_map[line];
		//add to index
		index.add(1, sample_vector);
	}
	cout << "...added " << index.ntotal << " samples to the index." << endl;
	db_file_stream.close();
	return index;
}


faiss::IndexHNSWFlat build_hnsw_index(string database_IDs, string database_data, char delim, int num_elements){
	// make map of sample ID to sample encoding (embedding)
	map<string, float*> ID_database_map = make_ID_data_map(database_data, delim, num_elements);

	faiss::IndexHNSWFlat index(num_elements, 40);
	index.hnsw.efSearch = 40;
	
	//int i = 0;
	//float* sample_arr = new float[num_elements];

	// read through databse file (list of samples in database)
	// and pull those vectors from the map to add to index
	ifstream db_file_stream;
	db_file_stream.open(database_IDs);
	if (!db_file_stream.is_open()){
		cout << "Failed to open: " << database_IDs << endl;
		exit(1);
	}
	string line;
	while (getline(db_file_stream, line)){
		// pull sample vector from map
		float* sample_vector = ID_database_map[line];
		//add to index
		index.add(1, sample_vector);
	}
	cout << "...added " << index.ntotal << " samples to the index." << endl;
	db_file_stream.close();
	return index;
}


map<string, float*> make_ID_data_map(string data_file, char delim, int num_elements){
	/*
	 * Returns a map where key is sampleID and 
	 * value is a vector of floats for
	 * the encoding (or embedding)
	 */

	// open data_file
	ifstream data_file_stream;
	data_file_stream.open(data_file);
	if( !data_file_stream.is_open()){
		cout << "Failed to open: " << data_file << endl;
		exit(1);
	}

	// returned map
	map<string, float*> ID_data_map;

	// read file
	string line;
	while (getline(data_file_stream, line)){
		size_t start;
		size_t end = 0;
		
		//vector<float> sample_vector;
		string sampleID = "";
		float* sample_vector = new float[num_elements];
		int float_i = 0;
		while ((start = line.find_first_not_of(delim, end)) != std::string::npos){
			end = line.find(delim, start);
			if (sampleID == ""){
				sampleID = line.substr(start, end - start);
			}else{
				float f = stof(line.substr(start, end - start));
				sample_vector[float_i] = f;
				float_i++;
			}
		}
		pair<string, float*> p (sampleID, sample_vector);
   	        ID_data_map.insert(p);
	}	

	return ID_data_map;

}

