#include <iostream>
#include <iterator>
#include <sstream>
#include <string>
#include <vector>

#include "utils.h"

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


// split a string into a vector
//vector<float> split(){
//	vector<float> vecFloats;
//
//}
//
//string s, tmp; 
//stringstream ss(s);
//vector<string> words;
//
//while(getline(ss, tmp, ',')){
//    words.push_back(tmp);
//    .....
//}
