#ifndef FAISSPM_H
#define FAISSPM_H

#endif //FAISSPM_H

#include <iostream>
#include <string>
#include <faiss/IndexFlat.h>

using namespace std;
using idx_t = faiss::Index::idx_t;

int faissMain(string eFile, int numV, int numS, int numQ);
