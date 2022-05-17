#include <algorithm>
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <faiss/IndexFlat.h>

#include "faiss_pm.h"

using namespace std;
/**
 * Copyright (c) Facebook, Inc. and its affiliates.
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


// 64-bit int
using idx_t = faiss::Index::idx_t;


int ss(float* database, float* queries, int numSamples, int numVariants, int numQueries){ 
	//int d = numVariants;      // dimension
	//int nb = numSamples; // database size
	//int nq = numQueries;  // nb of queries


	faiss::IndexFlatL2 index(numVariants); // call constructor
	printf("is_trained = %s\n", index.is_trained ? "true" : "false");
	index.add(numSamples, database); // add vectors to the index
	printf("ntotal = %zd\n", index.ntotal);

//	for (int i = 0; i < (d*nb); i++){
//		cout << database[i];
//	}
	
	int k = 4; // number of nearest neightbors to return
	
	// FAISS on database with queries
	{
		idx_t* I = new idx_t[k * numQueries];
		float* D = new float[k * numQueries];

		index.search(numQueries, queries, k, D, I);

		// print results
		cout << "RESULTS FOR SS:" << endl;
		for (int i = 0; i < numQueries; i++){
			cout << "  Query " << i << ": ";
			for (int j = 0; j < k; j++){
				cout << I[i * k * j] << "\t";
			}
			cout << endl;
		}
		
		delete[] I;
		delete[] D;
	}
	
	delete[] database;
	delete[] queries;

	return 0;
}

