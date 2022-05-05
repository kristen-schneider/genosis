#include <iostream>
#include <string>
#include <vector>
#include "utils.h"

using namespace std;

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
