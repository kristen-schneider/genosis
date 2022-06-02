#include <fstream>
#include <iostream>
#include <math.h>
#include <sstream>
#include <string>

#include "metrics.h"
#include "readEncoding.h"

using namespace std;

/*
 * compute euclidean distance
 * sqrt (sumN(v1-v2)^2)
 */
float euclidean_distance(float* vec1, float* vec2, int segLength){

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

/*
 * counts number of mismatches betweeen two vectors 
 */
float exact_match(float* vec1, float* vec2, int segLength){

        float numMismatches = 0;
        for (int i = 0; i < segLength; i++){
                if(vec1[i] != vec2[i]){
                        numMismatches++;
                }
        }
        return numMismatches;
}

/*
 * counts number of shared nonreference genotypes
 * so if they are both het or het and hom alt += one
 * the count of alleles where they share a non-refrence allele
 */
float sharedNRG(float* vec1, float* vec2, int segLength){
	float numNRG = 0;
	
	for (int i = 0; i < segLength; i++){
                if(vec1[i] != 0 && vec2[i] != 0){
                        numNRG++;
                }
        }

	return numNRG;
}
