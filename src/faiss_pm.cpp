#include <algorithm>
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>
//#include <math.h>
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


int ss(float* xb, float* xq, int numSamples, int numVariants, int numQueries){ 
	int d = numVariants;      // dimension
	int nb = numSamples; // database size
	int nq = numQueries;  // nb of queries

	//std::mt19937 rng;
	//std::uniform_real_distribution<> distrib;

	//float* xb = new float[d * nb];
	//float* xq = new float[d * nq];

	//float xb[nb][d] = {{1.f, 1.f, 1.f, 1.f, 1.f}, {2, 2, 2, 2, 2}, {3, 3, 3, 3, 3}};

	//for (int i = 0; i < nb; i++) {
	//	for (int j = 0; j < d; j++)
	//		xb[d*i+j] = 1.;
	//}

	//float xq[numVariants * nq] = {0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 1};

	faiss::IndexFlatL2 index(d); // call constructor
	printf("is_trained = %s\n", index.is_trained ? "true" : "false");
	index.add(nb, xb); // add vectors to the index
	printf("ntotal = %zd\n", index.ntotal);

//	for (int i = 0; i < (d*nb); i++){
//		cout << xb[i];
//	}
	
	int k = 4;
	
	{ // sanity check: search 5 first vectors of xb
		idx_t* I = new idx_t[k * nq];
		float* D = new float[k * nq];
		index.search(nq, xb, k, D, I);

	//	// print results for sanity check
    	//	cout << "RESULTS 1:" << endl;
	//	for (int i = 0; i < nq; i++){
	//		for (int j = 0; j < k; j++){
	//			cout << I[i * k * j];
	//		}
	//		cout << endl;
	//	}
			
        	delete[] I;
		delete[] D;
	}

	{
		idx_t* I = new idx_t[k * nq];
		float* D = new float[k * nq];

		index.search(nq, xq, k, D, I);

		// print results
		cout << "RESULTS FOR SS:" << endl;
		for (int i = 0; i < nq; i++){
			cout << "  Query " << i << ": ";
			for (int j = 0; j < k; j++){
				cout << I[i * k * j] << " ";
			}
			cout << endl;
		}
		
		delete[] I;
		delete[] D;
	}
	
	delete[] xb;
	delete[] xq;

	return 0;
}
