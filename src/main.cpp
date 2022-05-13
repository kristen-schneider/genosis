#include <iostream>

//#include "readVCF.h"
#include "readEncoding.h"
#include "faiss_pm.h"

using namespace std;

// code to read VCF and write to encoded file is commented out
// only code for reading encoded file and performing FAISS will be run
int main(void){

	// temp variables
	int numSamples = 3;
	int numVariants = 9;
	int numQueries = 3;

	string encodingtxt = "/home/sdp/precision-medicine/data/encoded/small_encoding.txt";

	//cout << "Start of encoding." << endl;

	//cout << "Reading VCF file." << endl;
	//sliceVCF();

	cout << "Reading Encoded file." << endl;
	float* xb = read_test(encodingtxt, numSamples, numVariants);
	cout << "Done Reading Encoded file." << endl;
	
	// queries
	float seed[numVariants * numQueries] = {0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 1.f, 1.f, 1.f, 1.f, 1.f, 1.f, 1.f, 1.f, 1.f, 1.f, 2.f, 2.f, 2.f, 2.f, 2.f, 2.f, 2.f, 2.f, 2.f, 2.f}; 
	float* xq = new float[numVariants * numQueries];
	for (int i = 0; i < (numVariants * numQueries); i++) {
        	xq[i] = seed[i];
	}


	cout << endl << "Starting FAISS." << endl;
	ss(xb, xq, numSamples, numVariants, numQueries);
	cout << "End of FAISS." << endl;

	return 0;
}
