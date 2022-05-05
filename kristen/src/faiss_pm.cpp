#include <algorithm>
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>
#include <vector>
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

float* faiss_flat(int numSamples, int numVariants, string* cohort_arr, float* cohort_float_arr){
	
	//float *cohort_float_arr[numSamples];
	
	for (int i = 0; i < numSamples; i++) {
		float sample_float_arr[numVariants];
		// convert a string to an array of floats
		istringstream ss(cohort_arr[i] );
  		copy(
    			istream_iterator <float> ( ss ),
    			istream_iterator <float> (),
    			sample_float_arr
    		);

		cohort_float_arr[i] = sample_float_arr;
		cout << cohort_arr[i] << endl;
	}
	return cohort_float_arr;	
}


int faiss_tutorial() { 
    int d = 64;      // dimension
    int nb = 100000; // database size
    int nq = 10000;  // nb of queries

    std::mt19937 rng;
    std::uniform_real_distribution<> distrib;

    float* xb = new float[d * nb];
    float* xq = new float[d * nq];

    for (int i = 0; i < nb; i++) {
        for (int j = 0; j < d; j++)
            xb[d * i + j] = distrib(rng);
        xb[d * i] += i / 1000.;
    }

    for (int i = 0; i < nq; i++) {
        for (int j = 0; j < d; j++)
            xq[d * i + j] = distrib(rng);
        xq[d * i] += i / 1000.;
    }

    faiss::IndexFlatL2 index(d); // call constructor
    printf("is_trained = %s\n", index.is_trained ? "true" : "false");
    index.add(nb, xb); // add vectors to the index
    printf("ntotal = %zd\n", index.ntotal);

    int k = 4;

    { // sanity check: search 5 first vectors of xb
        idx_t* I = new idx_t[k * 5];
        float* D = new float[k * 5];

        index.search(5, xb, k, D, I);

        // print results
        printf("I=\n");
        for (int i = 0; i < 5; i++) {
            for (int j = 0; j < k; j++)
                printf("%5zd ", I[i * k + j]);
            printf("\n");
        }

        printf("D=\n");
        for (int i = 0; i < 5; i++) {
            for (int j = 0; j < k; j++)
                printf("%7g ", D[i * k + j]);
            printf("\n");
        }

        delete[] I;
        delete[] D;
    }

    { // search xq
        idx_t* I = new idx_t[k * nq];
        float* D = new float[k * nq];

        index.search(nq, xq, k, D, I);

        // print results
        printf("I (5 first results)=\n");
        for (int i = 0; i < 5; i++) {
            for (int j = 0; j < k; j++)
                printf("%5zd ", I[i * k + j]);
            printf("\n");
        }

        printf("I (5 last results)=\n");
        for (int i = nq - 5; i < nq; i++) {
            for (int j = 0; j < k; j++)
                printf("%5zd ", I[i * k + j]);
            printf("\n");
        }

        delete[] I;
        delete[] D;
    }

    delete[] xb;
    delete[] xq;

    return 0;
}
