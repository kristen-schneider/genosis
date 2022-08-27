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

	cout << "Running FAISS on one segment encoding file." << endl;
	
	string encoding_file = argv[1];
	string query_file = argv[2];
	int k = stoi(argv[3]);
	// convert delim to char for differnt delims
	string s_delim = argv[4];
	char delim = '\0';
	if (s_delim == "space"){ delim = ' ';}


	int encoding_dimensions[2];	// [num_samples, num_variants]
	int query_dimensions[2];	// [num_queries, num_variants]
	get_dimensions(encoding_file, encoding_dimensions, delim);
	get_dimensions(query_file, query_dimensions, delim);
	int num_variants = encoding_dimensions[0];
	int num_samples = encoding_dimensions[1];
	int num_queries = query_dimensions[1];

	cout << "Number samples: " << num_samples << endl;
	cout << "Number variants: " << num_variants << endl;
	cout << "Number queries: " << num_queries << endl;
	


//	cout << "Starting FAISS for " << encoding_file << endl;
//
//	auto start = high_resolution_clock::now();
//	cout << "Building index..." << endl;
//	/*
//	faiss::IndexHNSWFlat index = build_faiss_index(encoding_file, \
//                        num_variants, \
//                        num_samples);
//	*/
//	faiss::IndexFlatL2 index = build_faiss_index(encoding_file, \
//			num_variants, \
//			num_samples, 
//			delim);
//
//	cout << "Running search..." << endl;
//	search(index, k, query_file, \
//		 	num_queries, \
//			num_variants, \
//			delim);	
//		
//	auto stop = high_resolution_clock::now();
//	auto duration_file = duration_cast<microseconds> (stop - start);
//	cout << "TIME:file:" << duration_file.count() << endl;
}
