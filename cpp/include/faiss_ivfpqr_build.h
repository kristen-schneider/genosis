#ifndef FAISS_IVFPQR_H
#define FAISS_IVFPQR_H

#endif //FAISS_IVFPQR_H

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

faiss::IndexIVF build_write_ivfpqr_index(string database_IDs, string database_encodings, const char* index_out_file);
//map<string, float*> make_ID_data_map(string data_file, char delim, int num_elements);
faiss::IndexIVF build_ivfpqr_index(int vector_size, string database_IDs, map<string, float*> ID_database_map);
