#include <iostream>

#include "readVCF.h"
#include "readEncoding.h"
#include "faiss_pm.h"

using namespace std;

//g++ main_encode.cpp slice.cpp utils.cpp -lhts -o encode

int main(void){

	// temp variables
	int numSamples = 10;
	int numVariants = 10;
	int numQueries = 3;

	cout << "Start of encoding." << endl;

	cout << "Reading VCF file." << endl;
	//sliceVCF();

	cout << "Reading Encoded file." << endl;
	//float* xb = read_test(numSamples, numVariants);
	
	float seed[numVariants * numQueries] = {0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 1.f, 1.f, 1.f, 1.f, 1.f, 1.f, 1.f, 1.f, 1.f, 1.f, 2.f, 2.f, 2.f, 2.f, 2.f, 2.f, 2.f, 2.f, 2.f, 2.f}; 
	float* xq = new float[numVariants * numQueries];
	for (int i = 0; i < (numVariants * numQueries); i++) {
        	xq[i] = seed[i];
	}

	cout << "FAISS." << endl;
	//ss(xb, xq, numSamples, numVariants, numQueries);

	cout << "End of encoding." << endl;
	return 0;
}
