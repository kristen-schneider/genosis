#ifndef NONREFGENOTYPES_H
#define NONREFGENOTYPES_H

#endif //NONREFGENOTYPES_H

#include <iostream>
#include <map>
#include <string> 

using namespace std;

map<int, int> count_non_reference_genotypes(string queriesTxt, string encodingTxt, int nQ, int nS, int nV);
map<int, int> encoding_count(string encodingTxt, float* qArr, int nS, int nV);
int count_nonref_genotypes(int sample, int nV, float* S, float* Q);
