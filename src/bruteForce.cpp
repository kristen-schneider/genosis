#include <iostream>
#include <math.h>

#include "bruteForce.h"

using namespace std;

float euclidean_distance(const float* vec1, float* vec2, int segLength){

	for (int i = 0; i < segLength; i++)
                cout << vec1[i] << " ";
	cout << endl;
	for (int i = 0; i < segLength; i++)
                cout << vec2[i] << " ";
	cout << endl;

	float eucDist = 0;
	float sum = 0;
	for (int i = 0; i < segLength; i++){
		float diff = vec1[0]-vec2[0];
		float diffSqrd = pow(diff, 2);
		sum += diffSqrd;
	}
	eucDist = sqrt(sum);

	return eucDist;
}

// returns an array from x[start:end]
float* arrSlice(float* x, int start, int end){
    float slice[end-start];
    int i_x = start;
    for (int i = 0; i < (end-start); i++){
        slice[i] = x[i_x];
        i_x ++;
    }
    return slice;
}
