#ifndef FAISS_L2_SEARCH
#define FAISS_L2_SEARCH

#endif //FAISS_L2_SEARCH

#include <cstdlib>
#include <faiss/IndexFlat.h>
#include <faiss/IndexHNSW.h>
#include <faiss/index_io.h>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>


// 64-bit int
using idx_t = faiss::Index::idx_t;
using namespace std;


float* make_queries_arr(string query_data, string query_IDs, char delim, int num_queries, int num_elements);
vector<string> make_query_ID_vector(string query_IDs);
