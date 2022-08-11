#include <iostream>
#include <fstream>
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

// return number of samples and number of variants in encoding file
void get_dimensions(string encodedTXT, int* dimensions){
	
	ifstream eFile;
	eFile.open(encodedTXT);
	if (!eFile.is_open()){
		cout << "Failed to open: " << encodedTXT << endl;
	}else{
		int num_samples = 0;
		int num_variants = 0;
		string line;
		while (getline(eFile, line)){
			if (num_variants == 0){
				num_variants = line.length();;
			}
			num_samples ++;	
		}
		dimensions[0] = num_samples;
		dimensions[1] = num_variants;
	}
}


