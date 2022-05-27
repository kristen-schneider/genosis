#include <algorithm>
#include <faiss/IndexFlat.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include "faissPM.h"


using namespace std;
/**
 * Copyright (c) Facebook, Inc. and its affiliates.
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


// 64-bit int
using idx_t = faiss::Index::idx_t;

int faissMain(string encodedFile){
	
	// ifstream to encoded file
        ifstream inFile;
	// open encoded file
        inFile.open(encodedFile);
        if ( !inFile.is_open() ) {
                cout << "Failed to open: " << encodedFile << endl;
        }

	// read encoded file line by line
	string line;
	if(inFile.is_open()){
                while(getline(inFile, line)){
			int segLength = line.length();
			string s;
			float f;

			// convert string line to float array
			float* lineArr = new float[segLength];
			for (int c = 0; c < segLength; c++){
				s = line[c];
				f = stof(s);
				lineArr[c] = f;	
			}

			// add array to index

			
		}
	}
	// closed encoded file
	inFile.close();
	return 0;

}
