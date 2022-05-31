#ifndef BUILDINDEX_H
#define BUILDINDEX_H

#endif //BUILDINDEX_H

#include <iostream>
#include <string>
#include <faiss/IndexFlat.h>

using namespace std;
using idx_t = faiss::Index::idx_t;

faiss::IndexFlatL2 build_faiss_index(string eFile, int numV, int numS);
faiss::IndexFlatL2 build_faiss_index_segments(string eFile, int start, int lengthSeg, int numS);
