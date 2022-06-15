#include <faiss/IndexFlat.h>
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

faiss::IndexFlatL2 build_faiss_index_segments(string encodedFile, int start, int lengthSegment, int numSamples){
	cout << "INDEX_FLAT_L2" << endl;
	// setup for FAISS
	faiss::IndexFlatL2 index(lengthSegment);
	if (index.is_trained == 1){cout << "...index is trained." << endl;}
	else{cerr << "...INDEX IS NOT TRAINED." << endl;}
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
			index.add(1, singleVector);	
			delete[] singleVector;
			lineCount++;
		}

	}
	cout << "...added " << index.ntotal << " vectors to index." << endl;
	// closed encoded file
	inFile.close();
	inFile.seekg(0);
	inFile.clear();
	return index;
}

faiss::IndexFlatIP build_faiss_index_segments_IP(string encodedFile, int start, int lengthSegment, int numSamples){
	cout << "INDEX_FLAT_IP" << endl;
	// setup for FAISS
        faiss::IndexFlatIP index(lengthSegment);
        if (index.is_trained == 1){cout << "...index is trained." << endl;}
        else{cerr << "...INDEX IS NOT TRAINED." << endl;}
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
                        index.add(1, singleVector);
                        delete[] singleVector;
                        lineCount++;
                }

        }
	cout << "...added " << index.ntotal << " vectors to index." << endl;
        // closed encoded file
        inFile.close();
        inFile.seekg(0);
        inFile.clear();
        return index;
}




