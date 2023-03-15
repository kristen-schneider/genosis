#include <cstdlib>
#include <faiss/IndexFlat.h>
#include <faiss/IndexHNSW.h>
#include <faiss/IndexPQ.h>
#include <faiss/IndexScalarQuantizer.h>
#include <faiss/AutoTune.h>
#include <faiss/index_factory.h>
#include <faiss/index_io.h>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>

#include "faiss_ivfpqr_build.h"
#include "faiss_utils.h"


// 64-bit int
using idx_t = faiss::Index::idx_t;
using namespace std;


int main(int argc, char* argv[]){
	
	string database_IDs = argv[1];
	string databse_encodings = argv[2];
	const char* index_out_file = argv[3];


	// FAISS IVFPQR index parameters
	int nlist = 100;
	int nprobe = 10;
	int m = 8;
	int ksub = 8;
	int nbits = 8;
	int niter = 10;
	int nredo = 5;
	int verbose = 2;
	
	// Build FAISS IVFPQR index
	faiss::Index *index = faiss_ivfpqr_build(database_IDs, databse_encodings, nlist, nprobe, m, ksub, nbits, niter, nredo, verbose);
	// write index to file
	faiss::write_index(index, index_out_file);
	delete index;
	return 0;

	//faiss::IndexIVF faiss_ivfpqr_index;
	//faiss_ivfpqr_index = build_write_ivfpqr_index(database_IDs, databse_encodings, index_out_file);
}

