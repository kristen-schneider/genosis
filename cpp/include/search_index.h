#ifndef SEARCHINDEX_H
#define SEARCHINDEX_H

#endif //SEARCHINDEX_H

#include <iostream>
#include <string>
#include <faiss/IndexFlat.h>
#include <faiss/IndexHNSW.h>

using namespace std;
using idx_t = faiss::Index::idx_t;

void search(const faiss::IndexFlatL2 &index, int k, string database_IDs, string query_IDs, string query_data, int num_queries, int num_variants, char delim);
//void search(const faiss::IndexHNSWFlat &index, int k, string queriesTXT, int num_queries, int num_variants);
float* make_queries_arr(string query_data, string query_IDs, char delim, int num_queries, int num_elements);
vector<string> make_query_ID_vector(string query_IDs);
