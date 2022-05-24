#include <iostream>
#include <math.h>

#include "bruteForce.h"

using namespace std;

float euclidean_distance(float* vec1, float* vec2, int segLength){

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
