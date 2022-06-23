#ifndef SEARCHINDEX_H
#define SEARCHINDEX_H

#endif //SEARCHINDEX_H

#include <iostream>
#include <string>
#include <faiss/IndexFlat.h>
#include <faiss/IndexHNSW.h>

using namespace std;
using idx_t = faiss::Index::idx_t;

void similarity_search(faiss::IndexFlatL2 index, string qFile, int start, int segLength, int numV, int numS, int numQ, int k, string txtName);
//void similarity_search(faiss::IndexHNSWFlat index, string qFile, int start, int segLength, int numV, int numS, int numQ, int k, string txtName);
