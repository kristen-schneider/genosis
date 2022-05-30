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

	// path to encoded file
	string encodingtxt = "/home/sdp/precision-medicine/data/encoded/new.encoded.txt";//ALL.wgs.svs.genotypes.encoded.txt";
	string queriestxt = "/home/sdp/precision-medicine/data/queries/new.queries.txt";//ALL.wgs.svs.genotypes.queries.txt";

	// DONE. Start FAISS..
	cout << "---------" << endl;
	cout << "Starting similarity searching using FAISS..." << endl;
	
	cout << "\n1.Builing index for " << encodingtxt << "..." << endl;
	faiss::IndexFlatL2 index = build_faiss_index(encodingtxt, numVariants, numSamples, numQueries);
	cout << "\n2.Running similairty search..." << endl;
	similarity_search(index, queriestxt, numVariants, numSamples, numQueries, k);	
	
	cout << "End of FAISS." << endl;
	cout << "---------" << endl;
	cout << "Starting Brute Force." << endl;
	int x = brute_force_main(encodingtxt, queriestxt, numVariants, numSamples, numQueries);
	

	return 0;
}
