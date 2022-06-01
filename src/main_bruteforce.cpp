#include <iostream>
#include <faiss/IndexFlat.h>

#include "buildIndex.h"
#include "searchIndex.h"
#include "bruteForce.h"

using namespace std;
// 64-bit int
using idx_t = faiss::Index::idx_t;

// code to read VCF and write to encoded file is commented out
// only code for reading encoded file and performing FAISS will be run
int main(void){

	// MAKE CHANGES TO THESE VARIABLES 
	// ...to be automated later...
	
	int numVariants = 2548903;//68819; // number of variants (cols) in encoding.txt
	int numSamples = 2548; // number of samples (rows) in encoding.txt
	int numQueries = 1; // number of queries
	int k = 2548;
	int segmentLength = 500; // length of a single vector
	// path to encoded file
	string encodingtxt = "/home/sdp/precision-medicine/data/encoded/new.encoded.txt";//ALL.wgs.svs.genotypes.encoded.txt";
	string queriestxt = "/home/sdp/precision-medicine/data/queries/new.queries.txt";//ALL.wgs.svs.genotypes.queries.txt";
	
	/*
	int numVariants = 9;
	int numSamples = 15;
	int numQueries = 1; // number of queries
	int k = 15;
	int segmentLength = 3; // length of a single vector

	// path to encoded file
	string encodingtxt = "/home/sdp/precision-medicine/data/encoded/test.encoded.txt";//ALL.wgs.svs.genotypes.encoded.txt";
	string queriestxt = "/home/sdp/precision-medicine/data/queries/test.queries.txt";//ALL.wgs.svs.genotypes.queries.txt";
	*/	
	
	int numSegments = numVariants/segmentLength;// + (numVariants % segmentLength != 0);
	cout << numSegments << endl;


	cout << "Starting Brute Force." << endl;
	int start_bf = 0;
	for (int i = 0; i < numSegments; i ++){
		cout << "\nSegment: " << start_bf << "-" << start_bf+segmentLength << endl;
		int x = brute_force_main(encodingtxt, queriestxt, start_bf, segmentLength, numVariants, numSamples, numQueries, numSegments);
		start_bf += segmentLength;
	}
	
	if (numVariants % segmentLength != 0){
		int lastSegmentLength = numVariants - (numSegments * segmentLength);
		cout << "\nLAST SEG DIFF";
		cout << "\nSegment: " << start_bf << "-" << start_bf+lastSegmentLength << endl;
		int x = brute_force_main(encodingtxt, queriestxt, start_bf, lastSegmentLength, numVariants, numSamples, numQueries, numSegments);	
	}

	

	return 0;
}
