#ifndef METRICS_H
#define METRICS_H

#endif //METRICS_H

#include <iostream>
#include <string>

using namespace std;

//int brute_force_main(string eFile, string qFile, int start, int lenQuery, int numV, int numS, int numQ, int numSeg);
//float* compute_one_query(float* q, string eFile, int start, int lenQUery, int numV, int numS, int numQ);
float euclidean_distance(float* vec1, float* vec2, int segLength);
float exact_match(float* vec1, float* vec2, int segLength);
float sharedNRG(float* vec1, float* vec2, int segLength);
