#include <iostream>

#include "readEncoding.h"
#include "nonR

using namespace std;

int count_non_reference_genotypes(string queriesTxt, string encodingTxt, int nQ, nV, nS){

	nonReferenceCount = 0;

	// path to queries file
        ifstream qFile;

        // open queries data (.txt file)
        qFile.open(queriesTxt);
        if ( !qFile.is_open() ) {
                cout << "Failed to open: " << queriesTxt << endl;
        }

        // to store all queries
        float* queriesArr = new float[nV * nQ];
        int Q = 0;
        if(qFile.is_open()){
                string line;
                while(getline(qFile, line)){

	return nonReferenceCount;
}
