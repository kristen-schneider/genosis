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

#include "faiss_l2_build.h"
#include "faiss_utils.h"


// 64-bit int
using idx_t = faiss::Index::idx_t;
using namespace std;


int main(int argc, char* argv[]){
	
	string database_IDs = argv[1];
	string databse_encodings = argv[2];
	const char* index_out_file = argv[3];

	faiss::IndexFlatL2 faiss_l2_index;
	faiss_l2_index = build_write_l2_index(database_IDs, databse_encodings, index_out_file);
	return 0;
}

faiss::IndexFlatL2 build_write_l2_index(string database_IDs, string database_encodings, const char* index_out_file){

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
	faiss::IndexFlatL2 faiss_l2_index = build_l2_index(vector_size, database_IDs, ID_encodings_map);
	cout << "writing index to " << index_out_file << endl;
	faiss::write_index(&faiss_l2_index, index_out_file);
	return faiss_l2_index;
}


faiss::IndexFlatL2 build_l2_index(int vector_size, string database_IDs, map<string, float*> ID_encodings_map){
//const FaissIndex* build_l2_index(int vector_size, string database_IDs, map<string, float*> ID_encodings_map, string index_file){
	cout << "building faiss l2 index." << endl;
	
	//const FaissIndex* index(vector_size);
	faiss::IndexFlatL2 index(vector_size);
	
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
	
	cout << "mapping encoding vectors to samples" << endl;
	
	// read file
	string line;
	while (getline(data_file_stream, line)){
		size_t start;
		size_t end = 0;

		//vector<float> sample_vector;
		string sampleID;
		auto* sample_vector = new float[num_elements];
		int float_i = 0;
		while ((start = line.find_first_not_of(delim, end)) != std::string::npos){
			end = line.find(delim, start);
			if (sampleID.empty()){
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
	cout << "finished mapping encoding vectors" << endl;
	return ID_data_map;
}
