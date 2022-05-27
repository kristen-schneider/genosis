#include <algorithm>
#include <faiss/IndexFlat.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include "faissPM.h"


/**
 * Copyright (c) Facebook, Inc. and its affiliates.
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


// 64-bit int
using idx_t = faiss::Index::idx_t;
using namespace std;

int faissMain(string encodedFile, int numVariants, int numSamples, int numQueries){

	// setup for FAISS
	faiss::IndexFlatL2 index(numVariants);
	printf("is_trained = %s\n", index.is_trained ? "true" : "false");

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
			int segLength = line.length();
			if (segLength != numVariants){
				cerr << "\t!! segment length does not equal number of variants !!" << endl;
			}
			string s;
			float f;

			// convert string line to float array
			float* singleVector = new float[segLength];
			for (int c = 0; c < segLength; c++){
				s = line[c];
				f = stof(s);
				singleVector[c] = f;	
			}
			// add array to index
			index.add(1, singleVector);	
			cout << "added vector: " << lineCount << endl;
			lineCount++;
			delete[] singleVector;
		}

	}
	printf("ntotal = %zd\n", index.ntotal);
	// closed encoded file
	inFile.close();
	return 0;

}
