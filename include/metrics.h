#ifndef METRICS_H
#define METRICS_H

#endif //METRICS_H

#include <iostream>
#include <string>

using namespace std;

float euclidean_distance(float* vec1, float* vec2, int segLength);
float mismatch(float* vec1, float* vec2, int segLength);
float smart_mismatch(float* vec1, float* vec2, int segLength);
float sharedNRG(float* vec1, float* vec2, int segLength);
float sharedNRGWeighted(float* vec1, float* vec2, int segLength);
