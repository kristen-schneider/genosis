#include <faiss/IndexFlat.h>
#include <fstream>
#include <iostream>
#include <math.h>
#include <sstream>
#include <string>

#include "searchIndex.h"
#include "readEncoding.h"

/**
 * Copyright (c) Facebook, Inc. and its affiliates.
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


// 64-bit int
using idx_t = faiss::Index::idx_t;
using namespace std;


void similarity_search(faiss::IndexFlatL2 index, string queriesFile, int lengthQuery, int numSamples, int numQueries, int k, string txtName){
	
	idx_t* I = new idx_t[k * numQueries];
	float* D = new float[k * numQueries];

	// get queries from file
	float* queries = read_queries(queriesFile, lengthQuery, numQueries);
	/*cout << "Query: " << endl;
	for (int i = 0; i < numVariants; i++){
        	cout << queries[i];
        }
	cout << endl;*/


	//  index.search(nq, xq, k, D, I);
	index.search(numQueries, queries, k, D, I);
	cout << "...search complete." << endl;

	// writing results
        cout << "...writing index results." << endl;
	ofstream outIndexFile;
	outIndexFile.open("/home/sdp/precision-medicine/data/txt/index."+txtName+".txt");
        for (int i = 0; i < numQueries; i++){
                for (int j = 0; j < k; j++){
                        outIndexFile << I[i * k + j] << "\t" << D[i * k + j] << endl;
                }
                outIndexFile << endl;
        }
	outIndexFile.close();
	
	/*
        cout << "..writing distance results." << endl;
	ofstream outDistanceFile;
	outDistanceFile.open("/home/sdp/precision-medicine/data/txt/distanceReults.txt");
        for (int i = 0; i < numQueries; i++){
                for (int j = 0; j < k; j++){
                        outDistanceFile << "\t" << sqrt(D[i * k + j]) << " ";
                }
                outDistanceFile << endl;
        }
	outDistanceFile.close();
	*/
	
	delete [] I;
	delete [] D;
}
