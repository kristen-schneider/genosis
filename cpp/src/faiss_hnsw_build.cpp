#include <cstdlib>
#include <faiss/IndexFlat.h>
#include <faiss/IndexHNSW.h>
#include <faiss/AutoTune.h>
#include <faiss/index_factory.h>
#include <faiss/index_io.h>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>

#include "faiss_hnsw_build.h"
#include "faiss_utils.h"


// 64-bit int
using idx_t = faiss::Index::idx_t;
using namespace std;


int main(int argc, char* argv[]){
	
	string database_IDs = argv[1];
	string databse_encodings = argv[2];
	const char* index_out_file = argv[3];

	faiss::IndexHNSWFlat faiss_hnsw_index;
	faiss_hnsw_index = build_write_hnsw_index(database_IDs, databse_encodings, index_out_file);
	return 0;
}

faiss::IndexHNSWFlat build_write_hnsw_index(string database_IDs, string database_encodings, const char* index_out_file){

	// delim
	char delim = ' ';

	// vector size
	int vector_size = count_length_input_vector(database_encodings, delim);
	cout << "size of vector: " << vector_size << endl;
	// number samples
	int number_samples = count_num_samples(database_IDs);
	cout << "number of samples: " << number_samples << endl;

	// make map of sample ID to sample encoding (embedding)
	map<string, float*> ID_encodings_map;
	ID_encodings_map = make_ID_data_map(database_encodings, delim, vector_size);

	// make faiss index
	faiss::IndexHNSWFlat faiss_hnsw_index = build_hnsw_index(vector_size, database_IDs, ID_encodings_map);
	cout << "writing index to " << index_out_file << endl;
	faiss::write_index(&faiss_hnsw_index, index_out_file);
	return faiss_hnsw_index;
}


faiss::IndexHNSWFlat build_hnsw_index(int vector_size, string database_IDs, map<string, float*> ID_encodings_map){
//const FaissIndex* build_hnsw_index(int vector_size, string database_IDs, map<string, float*> ID_encodings_map, string index_file){
	cout << "building faiss hnsw index." << endl;
	
	faiss::IndexHNSWFlat index(vector_size, 40);
	index.hnsw.efSearch = 40;

	// read through databse file (list of samples in database)
	// and pull those vectors from the map to add to index
	ifstream db_file_stream;
	db_file_stream.open(database_IDs);
	if (!db_file_stream.is_open()){
		cout << "Failed to open: " << database_IDs << endl;
		exit(1);
	}
	string line;
	cout << "...adding vectors" << endl;
	while (getline(db_file_stream, line)){
		// pull sample vector from map
		float* sample_vector = ID_encodings_map[line];
		//add to index
		index.add(1, sample_vector);
	}
	cout << "...added " << index.ntotal << " samples to the index." << endl;
	db_file_stream.close();
	
	// write index to file
	//const char* out_file = "out.out";
	//faiss::write_index(&index, index_file);
	return index;
}

