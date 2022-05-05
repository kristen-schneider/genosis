#include <iostream>
#include <math.h>

#include "faiss_ss"

void faiss(){
	int d = 64;                            // dimension
   	int nb = 100000;                       // database size
   	int nq = 10000;                        // nb of queries

	float *xb = new float[d * nb];
	float *xq = new float[d * nq];

	for(int i = 0; i < nb; i++) {
        	for(int j = 0; j < d; j++) xb[d * i + j] = drand48();
        	xb[d * i] += i / 1000.;
    	}
    	for(int i = 0; i < nq; i++) {
        	for(int j = 0; j < d; j++) xq[d * i + j] = drand48();
        	xq[d * i] += i / 1000.;
    	}
}
