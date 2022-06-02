#include <iostream>

#include "buildIndex.h"
#include "compare.h"
#include "metrics.h"
#include "searchIndex.h"

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
/*/
	
	
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
	int metric = -1;
/*
	cout << "------------------" << endl;
	cout << "Starting Euclidean Distance." << endl;
	int start_ed = 0;
	metric = 0;
	for (int i = 0; i < numSegments; i ++){
		cout << "\nSegment: " << start_ed << "-" << start_ed+segmentLength << endl;
		compare_main(encodingtxt, queriestxt, start_ed, segmentLength, numVariants, numSamples, numQueries, numSegments, metric);
		start_ed += segmentLength;
	}
	// last segment if differnt length	
	if (numVariants % segmentLength != 0){
		int lastSegmentLength = numVariants - (numSegments * segmentLength);
		cout << "\nSegment: " << start_ed << "-" << start_ed+lastSegmentLength << endl;
		compare_main(encodingtxt, queriestxt, start_ed, lastSegmentLength, numVariants, numSamples, numQueries, numSegments, metric);	
	}
	cout << "End of Euclidean Distance." << endl;
        cout << "------------------\n" << endl;
	

	cout << "------------------" << endl;
	cout << "Starting Counting Mismatches." << endl;
	int start_mm = 0;
	metric = 1;
	for (int i = 0; i < numSegments; i ++){
		cout << "\nSegment: " << start_mm << "-" << start_mm+segmentLength << endl;
		compare_main(encodingtxt, queriestxt, start_mm, segmentLength, numVariants, numSamples, numQueries, numSegments, metric);
		start_mm += segmentLength;
	}
	// last segment if differnt length	
	if (numVariants % segmentLength != 0){
		int lastSegmentLength = numVariants - (numSegments * segmentLength);
		cout << "\nSegment: " << start_mm << "-" << start_mm+lastSegmentLength << endl;
		compare_main(encodingtxt, queriestxt, start_mm, lastSegmentLength, numVariants, numSamples, numQueries, numSegments, metric);	
	}
	cout << "End of Counting Mismatches." << endl;
        cout << "------------------\n" << endl;
*/

	cout << "------------------" << endl;
	cout << "Starting Counting Non-reference Genotypes." << endl;
        int start_nrg = 0;
        metric = 2;
        for (int i = 0; i < numSegments; i ++){
                cout << "\nSegment: " << start_nrg << "-" << start_nrg+segmentLength << endl;
                compare_main(encodingtxt, queriestxt, start_nrg, segmentLength, numVariants, numSamples, numQueries, numSegments, metric);
                start_nrg += segmentLength;
        }
        // last segment if differnt length      
        if (numVariants % segmentLength != 0){
                int lastSegmentLength = numVariants - (numSegments * segmentLength);
                cout << "\nSegment: " << start_nrg << "-" << start_nrg+lastSegmentLength << endl;
                compare_main(encodingtxt, queriestxt, start_nrg, lastSegmentLength, numVariants, numSamples, numQueries, numSegments, metric);
        }
        cout << "End of Counting Mismatches." << endl;
        cout << "------------------\n" << endl;

	return 0;
}
