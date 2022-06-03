#include <iostream>
#include <faiss/IndexFlat.h>

#include "buildIndex.h"
#include "searchIndex.h"

using namespace std;
// 64-bit int
using idx_t = faiss::Index::idx_t;

// code to read VCF and write to encoded file is commented out
// only code for reading encoded file and performing FAISS will be run
int main(void){

	// MAKE CHANGES TO THESE VARIABLES 
	// ...to be automated later...
	/*
	int numVariants = 2548903;//68819; // number of variants (cols) in encoding.txt
	int numSamples = 2548; // number of samples (rows) in encoding.txt
	int numQueries = 1; // number of queries
	int k = 2548;
	int segmentLength = 500; // length of a single vector
	// path to encoded file
	//string encodingtxt = "/home/sdp/precision-medicine/data/encoded/new.encoded2.txt";//ALL.wgs.svs.genotypes.encoded.txt
	//string queriestxt = "/home/sdp/precision-medicine/data/queries/new.queries2.txt";//ALL.wgs.svs.genotypes.queries.txt
	*/
	
	int numVariants = 9;
	int numSamples = 3;
	int numQueries = 1; // number of queries
	int k = 3;
	int segmentLength = 3; // length of a single vector

	// path to encoded file
	string encodingtxt = "/home/sdp/precision-medicine/data/encoded/short.encoded.txt";//ALL.wgs.svs.genotypes.encoded.txt
	string queriestxt = "/home/sdp/precision-medicine/data/queries/short.queries.txt";//ALL.wgs.svs.genotypes.queries.txt
	
	
	int numSegments = numVariants/segmentLength;// + (numVariants % segmentLength != 0);
	cout << numSegments << endl;


	// DONE. Start FAISS..
	cout << "---------" << endl;
	cout << "Starting similarity searching using FAISS..." << endl;
	
	cout << "\nENCODED FILE: " << encodingtxt << "..." << endl;
	//faiss::IndexFlatL2 index = build_faiss_index(encodingtxt, numVariants, numSamples);
	
	int start = 0;
	for (int i = 0; i < numSegments; i ++){
		cout << "\nSegment: " << start << "-" << start+segmentLength << endl;
		cout << "-Building index." << endl;
		faiss::IndexFlatL2 s_index = build_faiss_index_segments(encodingtxt, start, segmentLength, numSamples);
		cout << "-Running similairty search." << endl;
		similarity_search(s_index, queriestxt, start, segmentLength, numVariants, numSamples, numQueries, k, to_string(start));
		start += segmentLength;
	}
	if (numVariants % segmentLength != 0){
		int lastSegmentLength = numVariants - (numSegments * segmentLength);
		cout << "\nLAST SEG DIFF";// << endl;
		cout << "\nSegment: " << start << "-" << start+lastSegmentLength << endl;
                cout << "-Building index." << endl;
                faiss::IndexFlatL2 s_index = build_faiss_index_segments(encodingtxt, start, lastSegmentLength, numSamples);
                cout << "-Running similairty search." << endl;
                similarity_search(s_index, queriestxt, start, lastSegmentLength, numVariants, numSamples, numQueries, k, to_string(start));

	}
	/*
	cout << "\n2.Running similairty search..." << endl;
	similarity_search(index, queriestxt, numVariants, numSamples, numQueries, k);	
	*/
	
	cout << "End of FAISS." << endl;
	cout << "---------" << endl;


	return 0;
}
