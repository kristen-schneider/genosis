#include <chrono>
#include <cstdlib>
#include <iostream>
#include <faiss/IndexFlat.h>
#include <faiss/IndexHNSW.h>
#include <fstream>
#include <sstream>

#include "buildIndex.h"
#include "searchIndex.h"

using namespace std;
using namespace std::chrono;
// 64-bit int
using idx_t = faiss::Index::idx_t;

// code to read VCF and write to encoded file is commented out
// only code for reading encoded file and performing FAISS will be run
int main(int argc, char* argv[]){
	
	// path to encoded file
	string encodingtxt = argv[1];		// file with encoded data
	string queriestxt = argv[2];		// file with query data

	int numVariants = atoi(argv[3]);	// number of variants (cols)
	int numSamples = atoi(argv[4]);		// number of samples (rows)
	int numQueries= atoi(argv[5]); 		// number of queries sumbitted
	int k = atoi(argv[6]); 			// num nearest neighbors
	int segmentLength = atoi(argv[7]); 	// length of one segment


	int numSegments = numVariants/segmentLength;
	cout << "NUM SEGMENTS (floor): " << numSegments << endl;


	// DONE. Start FAISS..
	cout << "---------" << endl;
	cout << "Starting similarity searching using FAISS..." << endl;
	
	cout << "\nENCODED FILE: " << encodingtxt << "..." << endl;
	
	int start = 0;
	for (int i = 0; i < numSegments; i ++){
		auto startTime = high_resolution_clock::now();	
		cout << "\nSegment: " << start << "-" << start+segmentLength << endl;
		cout << "-Building index." << endl;
		
		faiss::IndexFlatL2 s_index = build_faiss_index_segments(encodingtxt, start, segmentLength, numSamples);
		//faiss::IndexHNSWFlat s_index = build_faiss_index_segments(encodingtxt, start, segmentLength, numSamples);
		cout << "-Running similairty search." << endl;
		similarity_search(s_index, queriestxt, start, segmentLength, numVariants, numSamples, numQueries, k, to_string(start));
		start += segmentLength;
		
		auto stopTime = high_resolution_clock::now();
		auto durationTime = duration_cast<microseconds>(stopTime - startTime);
		cout << "TIME: " << durationTime.count() << endl;
		cout << "end segment: " << i << endl;
	}
	if (numVariants % segmentLength != 0){
		int lastSegmentLength = numVariants - (numSegments * segmentLength);
		cout << "\nLAST SEG DIFF";// << endl;
		cout << "\nSegment: " << start << "-" << start+lastSegmentLength << endl;
                cout << "-Building index." << endl;

		faiss::IndexFlatL2 s_index = build_faiss_index_segments(encodingtxt, start, lastSegmentLength, numSamples);
		//faiss::IndexHNSWFlat s_index = build_faiss_index_segments(encodingtxt, start, lastSegmentLength, numSamples);
                cout << "-Running similairty search." << endl;
                similarity_search(s_index, queriestxt, start, lastSegmentLength, numVariants, numSamples, numQueries, k, to_string(start));

	}
	
	cout << "End of FAISS." << endl;
	cout << "---------" << endl;


	return 0;
}
