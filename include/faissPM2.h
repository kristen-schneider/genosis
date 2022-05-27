#ifndef FAISSPM2_H
#define FAISSPM2_H

#endif //FAISSPM2_H

#include <iostream>
#include <string>
#include <faiss/IndexFlat.h>

using namespace std;
using idx_t = faiss::Index::idx_t;

int similarity_search(faiss::IndexFlatL2 index, string qFile, int numV, int numS, int numQ);
