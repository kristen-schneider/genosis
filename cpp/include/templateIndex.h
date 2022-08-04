#ifndef TEMPLATE_INDEX_H
#define TEMPLATE_INDEX_H

#endif //TEMPLATE_INDEX_H

#include <iostream>
#include <string>

// 64-bit int
using idx_t = faiss::Index::idx_t;
using namespace std;

template <class indexType>
indexType buildIndex(indexType index, string encodedFile, int start, int lengthSegment, int numSamples);
