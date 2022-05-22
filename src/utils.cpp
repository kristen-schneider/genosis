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

/*
map<string, string> create_map(){
	map<string, string> m;
	m["key"] = "value";
	return m;
}
*/

/*
// initialize map for genotype encoding

map<string,int> create_map(){
  	map<string, int> m;
	m["0|0"] = 0;

	m["0|1"] = 1;
	m["1|0"] = 1;
	m["0|2"] = 1;
	m["2|0"] = 1;
	m["0|3"] = 1;
	m["3|0"] = 1;
	m["1|2"] = 1;
	m["2|1"] = 1;	
	m["1|3"] = 1;
	m["3|1"] = 1;

	m["1|1"] = 2;
	m["2|2"] = 2;
	m["3|3"] = 2;
	
	m[".|."] = 3;
	m["0|."] = 3;
	m[".|0"] = 3;
	m["1|."] = 3;
	m[".|1"] = 3;
	m["2|."] = 3;
	m[".|2"] = 3;
	m["3|."] = 3;
	m[".|3"] = 3;


	return m;
}
*/
