#include <cstdlib>
#include <faiss/IndexFlat.h>
#include <faiss/IndexHNSW.h>
#include <iostream>
#include <chrono>

#include "build_index.h"
#include "search_index.h"
#include "utils.h"

// 64-bit int
using idx_t = faiss::Index::idx_t;
using namespace std::chrono;
using namespace std;

int main(int argc, char* argv[]){

	string database_IDs = argv[1];
	string database_data = argv[2];
	string query_IDs = argv[3];
	string query_data = argv[4];
	int k = stoi(argv[5]);
	// convert delim to char for differnt delims
	string s_delim = argv[6];
	char delim = '\0';
	if (s_delim == "space"){ delim = ' ';}


	int num_elements = count_length_input_vector(database_data, delim);
	int num_db_samples = count_num_samples(database_IDs);
	int num_q_samples = count_num_samples(query_IDs);
	
	cout << "FAISS L2" << endl;
	cout << "Number database samples: " << num_db_samples << endl;
	cout << "Number query samples: " << num_q_samples << endl;
	cout << "Number vector elements: " << num_elements << endl;

	auto start = high_resolution_clock::now();
	cout << "Building Index..." << endl;
	faiss::IndexFlatL2 faiss_index = build_l2_index(database_IDs, database_data, delim, num_elements);
	cout << "Running search..." << endl;
	search_l2(faiss_index, k, database_IDs, query_IDs, query_data, num_q_samples, num_elements, delim);	
	auto stop = high_resolution_clock::now();
	auto duration_file = duration_cast<microseconds> (stop - start);
	cout << "TIME:file:" << duration_file.count() << endl;
}
