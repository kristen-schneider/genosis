#ifndef FAISS_PM_H
#define FAISS_PM_H

#endif //FAISS_PM_H

#include <iostream>
#include <string>
#include <faiss/IndexFlat.h>

using namespace std;
using idx_t = faiss::Index::idx_t;

int ss(float*, float*, int, int, int, int);
float FAISS_vs_BF(float* database, float* queries, int numSamples, int segmentLength, int numQueries, idx_t* I, float* D, int k);
