#ifndef COMPARE_H
#define COMPARE_H

#endif //COMPARE_H

#include <iostream>
#include <string>

using namespace std;

void compare_main(string eFile, string qFile, int start, int lenQuery, int numV, int numS, int numQ, int numSeg, int metric);
float* compute_one_query(float* q, string eFile, int start, int lenQUery, int numV, int numS, int numQ, int metric);
