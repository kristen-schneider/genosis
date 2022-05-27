#ifndef BRUTEFORCE_H
#define BRUTEFFORCE_H

#endif //BRUTEFFORCE_H

#include <iostream>
#include <string>

using namespace std;

int brute_force_main(string eFile, string qFile, int numV, int numS, int numQ);
float* compute_one_query(float* q, string eFile, int numV, int numS, int numQ);
float euclidean_distance(float* vec1, float* vec2, int segLength);
