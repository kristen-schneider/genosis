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

faiss::IndexFlatL2 build_index(string data_file, string database_file, char delim, int num_elements){
	// make map of sample ID to sample encoding (embedding)
	map<string, float*> ID_vector_map = make_ID_vector_map(data_file, delim, num_elements);
	
	faiss::IndexFlatL2 index(num_elements);
	//int i = 0;
	//float* sample_arr = new float[num_elements];

	// read through databse file (list of samples in database)
	// and pull those vectors from the map to add to index
	ifstream db_file_stream;
	db_file_stream.open(database_file);
	if (!db_file_stream.is_open()){
		cout << "Failed to open: " << database_file << endl;
		exit(1);
	}
	string line;
	while (getline(db_file_stream, line)){
		// pull sample vector from map
		float* sample_vector = ID_vector_map[line];
		//add to index
		index.add(1, sample_vector);
		delete[] sample_vector;
		//for (int j = 0; j < sample_vector.size(); j++){
			//float f = sample_vector.at(j);
			//sample_arr[i] = f;
			//i++;
		//}
	}
	cout << "...added " << index.ntotal << " samples to the index." << endl;
	db_file_stream.close();
	return index;
}

/*
 * Returns a map where key is sampleID and 
 * value is a vector of floats for
 * the encoding (or embedding)
 */
map<string, float*> make_ID_vector_map(string data_file, char delim, int num_elements){
	// open data_file
	ifstream data_file_stream;
	data_file_stream.open(data_file);
	if( !data_file_stream.is_open()){
		cout << "Failed to open: " << data_file << endl;
		exit(1);
	}

	map<string, float*> ID_vector_map;

	string line;
	size_t start;
	size_t end = 0;
	// read file
	while (getline(data_file_stream, line)){
		//vector<float> sample_vector;
		string sampleID = "";
		float* sample_vector = new float[num_elements];
		int float_i = 0;
		while ((start = line.find_first_not_of(delim, end)) != std::string::npos){
			if (sampleID == ""){
				sampleID = line.substr(start, end - start);
			}else{
				float f = stof(line.substr(start, end - start));
				sample_vector[float_i] = f;
				float_i++;
			}
		}
		pair<string, float*> p (sampleID, sample_vector);
   	        ID_vector_map.insert(p);
	}	

	return ID_vector_map;

}

//faiss::IndexHNSWFlat build_faiss_index(string encoding_file, int num_elements, int num_samples){
faiss::IndexFlatL2 build_faiss_index(string encoding_file, int num_elements, int num_samples, char delim){
	// setup for FAISS
	faiss::IndexFlatL2 index(num_elements);
	//faiss::IndexHNSWFlat index(num_elements, 64);

	// open encoded file
	ifstream eFile;
	eFile.open(encoding_file);
	if (!eFile.is_open()){
		cout << "Failed to open: " << encoding_file << endl;
		exit(1);
	}
	string line;
	while (getline(eFile, line)){
		string s;
		float f;

		// convert line to float array
		int i = 0;
		float* sample_vector = new float[num_elements];
		size_t start;
		size_t end = 0;
		if (delim == '\0'){
		 	for (int i = 0; i < num_elements; i ++){
				s = line[i];
				f = stof(s);
				sample_vector[i] = f;
			}
		}
		else{
			while ((start = line.find_first_not_of(delim, end)) != std::string::npos){
				end = line.find(delim, start);
				f = stof(line.substr(start, end - start));
				sample_vector[i] = f;
				i ++;
		
			}
		}
		/*
		for (int i = 0; i < num_elements; i ++){
			cout << sample_vector[i] << " ";
		}
		cout << endl;
		*/
		index.add(1, sample_vector);
		delete[] sample_vector;
	}
	cout << "...added " << index.ntotal << " samples to the index." << endl;
	eFile.seekg(0);
	eFile.close();
	eFile.clear();
	return index;
}

