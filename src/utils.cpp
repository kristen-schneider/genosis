#include <iostream>
//#include <iterator>
//#include <sstream>
//#include <fstream>
#include <string>
#include <vector>
#include <map>

//#include "utils.h"

using namespace std;


// transpose a vector of vectors of ints
vector<vector<int>> transpose(vector<vector<int>> &b){
	cout << "Transposing VMF to SMF..." << endl;
	// throw error if size is bad
	if (b.size() == 0) {
		cerr << "Error reaching to db." << endl;
	}
    	
	// transpose data
	vector<vector<int>> trans_vec(b[1].size(), vector<int>());
    	for (size_t i = 0; i < b.size(); i++) {
		for (size_t j = 0; j < b[i].size(); j++) {
            		trans_vec[j].push_back(b[i][j]);
        	}
    	}
    return trans_vec;
}

// returns an array from x[start:end]
float* arrSlice(float* x, int start, int end){
    float slice[end-start];
    int i_x = start;
    for (int i = 0; i < (end-start); i++){
        slice[i] = x[i_x];
        i_x ++;
    }
    return slice;
}
