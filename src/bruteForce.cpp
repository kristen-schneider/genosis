#include <iostream>
#include <math.h>

#include "bruteForce.h"

using namespace std;

float euclidean_distance(float* vec1, float* vec2, int segLength){

	for (int i = 0; i < segLength; i++)
                cout << vec1[i] << " ";
	cout << endl;
	for (int i = 0; i < segLength; i++)
                cout << vec2[i] << " ";
	cout << endl;

	float eucDist = 0;
	float sum = 0;
	for (int i = 0; i < segLength; i++){
		float diff = vec1[i]-vec2[i];
		float diffSqrd = pow(diff, 2);
		sum += diffSqrd;
	}
	eucDist = sqrt(sum);
	return eucDist;
}
