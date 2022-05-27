#ifndef READENCODINGS_H
#define READENCODINGS_H

#endif //READENCODINGS_H

#include <iostream>
#include <string> 

using namespace std;

float* read_encodings(string encodingtxt, int numSamples, int numVariants);
float* read_queries(string queriestxt, int numVariants, int numQueries);
