#include <iostream>
#include <map>

#include "nonRefGenotypes.h"

using namespace std;

int main(){
	
	string queriesTxt = "/home/sdp/precision-medicine/data/queries/test.queries.txt";
	string encodingTxt = "/home/sdp/precision-medicine/data/encoded/test.encoded.txt";
	int nQ = 1;
	int nS = 15;//2548;
	int nV = 9;//2548903;

	map<int, int> nonRefCount;
	nonRefCount = count_non_reference_genotypes(queriesTxt, encodingTxt, nQ, nS, nV);
	
	return 0;
}
