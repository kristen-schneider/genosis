#include <iostream>

#include "readVCF.h"
#include "readEncoding.h"
#include "faiss_pm.h"

using namespace std;

//g++ main_encode.cpp slice.cpp utils.cpp -lhts -o encode

int main(void){

	// temp variables
	int numSamples = 3;
	int numVariants = 9;

	cout << "Start of encoding." << endl;

	cout << "Reading VCF file." << endl;
	//sliceVCF();

	cout << "Reading Encoded file." << endl;
	read_test(numSamples, numVariants);

	cout << "FAISS." << endl;
	ss();

	cout << "End of encoding." << endl;
	return 0;
}
