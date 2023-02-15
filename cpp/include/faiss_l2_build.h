#ifndef FAISS_L2_H
#define FAISS_L2_H

#endif //FAISS_L2_H

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

faiss::IndexFlatL2 build_write_l2_index(string database_IDs, string database_encodings, const char* index_out_file);
//map<string, float*> make_ID_data_map(string data_file, char delim, int num_elements);
faiss::IndexFlatL2 build_l2_index(int vector_size, string database_IDs, map<string, float*> ID_database_map);
