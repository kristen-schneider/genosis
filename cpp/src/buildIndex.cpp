#include <faiss/IndexFlat.h>
#include <faiss/IndexHNSW.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include "buildIndex.h"


/**
 * Copyright (c) Facebook, Inc. and its affiliates.
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


// 64-bit int
using idx_t = faiss::Index::idx_t;
using namespace std;

faiss::IndexHNSWFlat build_faiss_index_segments(string encodedFile, int start, int lengthSegment, int numSamples){
//faiss::IndexFlatL2 build_faiss_index_segments(string encodedFile, int start, int lengthSegment, int numSamples){
	cout << "INDEX_HNSW_FLAT" << endl;
	
	// setup for FAISS
	faiss::IndexHNSWFlat IndexHNSWFlat(lengthSegment, 64);
	//faiss::IndexFlatL2 index(lengthSegment);
	
	//if (index.is_trained == 1){cout << "...index is trained." << endl;}
	//else{cerr << "...INDEX IS NOT TRAINED." << endl;}
	
	// ifstream to encoded file
        ifstream inFile;
	// open encoded file
        inFile.open(encodedFile);
        if ( !inFile.is_open() ) {
                cout << "Failed to open: " << encodedFile << endl;
        }

	// read encoded file line by line
	string line;
	int lineCount = 0;
	if(inFile.is_open()){
                while(getline(inFile, line)){
			string s;
			float f;
			// convert string line to float array
			float* singleVector = new float[lengthSegment];
			int i = 0;
			for (int c = start; c < start+lengthSegment; c++){
				s = line[c];
				f = stof(s);
				singleVector[i] = f;	
				i++;
			}
			/*cout << "adding vector: ";
			for (int i = 0; i < segLength; i++){
				cout << singleVector[i];
			}*/

			// add array to index
			IndexHNSWFlat.add(1, singleVector);	
			delete[] singleVector;
			lineCount++;
		}

	}
	cout << "...added " << IndexHNSWFlat.ntotal << " vectors to index." << endl;
	// closed encoded file
	inFile.seekg(0);
	inFile.close();
	inFile.clear();
	return IndexHNSWFlat;
}
