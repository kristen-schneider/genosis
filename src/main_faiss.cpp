#include <iostream>

//#include "readVCF.h"
#include "utils.h"
#include "readEncoding.h"
#include "chunkEncoding.h"
#include "faiss_pm.h"
#include "bruteForce.h"

using namespace std;

// code to read VCF and write to encoded file is commented out
// only code for reading encoded file and performing FAISS will be run
int main(void){

	// MAKE CHANGES TO THESE VARIABLES 
	// ...to be automated later...
	int numSamples = 5096;//2504; // number of samples (rows) in encoding.txt
	int numVariants = 2500;//68819; // number of variants (cols) in encoding.txt
	int numQueries = 5; // number of queries
	int segmentLength = 5;

	// path to encoded file
	string encodingtxt = "/home/sdp/precision-medicine/data/encoded/chr14.encoded.txt";//ALL.wgs.svs.genotypes.encoded.txt";
	string queriestxt = "/home/sdp/precision-medicine/data/queries/test.queries.txt";//ALL.wgs.svs.genotypes.queries.txt";

	// DONE. Start FAISS..
	// create an array which will holds queries
	float* xq = read_queries(queriestxt, numSamples, numVariants);

	cout << "Reading Encoded file." << endl;
	float* xb = read_encodings(encodingtxt, numSamples, numVariants);
	//float* chunk = split_encoding(encodingtxt, numSamples, segmentLength);
	cout << "Done Reading Encoded file." << endl;

	cout << endl << "Starting FAISS." << endl;
	ss(xb, xq, numSamples, numVariants, numQueries, segmentLength);
	cout << "End of FAISS." << endl;

	

	return 0;
}
