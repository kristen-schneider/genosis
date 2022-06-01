#ifndef EXACTMATCH_H
#define EXACTMATCH_H

#endif //EXACTMATCH_H

#include <iostream>
#include <string>

using namespace std;

int exact_match_main(string eFile, string qFile, int start, int lenQuery, int numV, int numS, int numQ, int numSeg);
float* compute_one_query(float* q, string eFile, int start, int lenQUery, int numV, int numS, int numQ);
float exact_match(float* vec1, float* vec2, int segLength);
