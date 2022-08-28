#ifndef BUILDINDEX_H
#define BUILDINDEX_H

#endif //BUILDINDEX_H

#include <iostream>
#include <faiss/IndexFlat.h>
#include <faiss/IndexHNSW.h>
#include <map>
#include <string>
#include <vector>

using namespace std;
using idx_t = faiss::Index::idx_t;

faiss::IndexFlatL2 build_index(string database_IDs, string database_data, char delim, int num_elements);
map<string, float*> make_ID_data_map(string data_file, char delim, int num_elements);


faiss::IndexFlatL2 build_faiss_index(string encodedTXT, int num_variants, int num_samples, char delim);
//faiss::IndexHNSWFlat build_faiss_index(string encodedTXT, int num_variants, int num_samples);

