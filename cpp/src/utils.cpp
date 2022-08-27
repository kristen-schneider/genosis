#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>

//#include "utils.h"

using namespace std;

// transpose a vector of vectors of ints
vector<vector<int>> transpose(vector<vector<int>> &vmf){
        /*
         * Takes a variant major format
         * vector of vector of ints
         * and transposes it to
         * sample major format.
         */
        cout << "Transposing VMF to SMF..." << endl;
        // throw error if size is bad
        if (vmf.size() == 0) {
                cerr << "Error reading variant major format." << endl;
        }

        // transpose data
        vector<vector<int>> smf_vec(vmf[0].size(), vector<int>());
        for (size_t i = 0; i < vmf.size(); i++) {
                for (size_t j = 0; j < vmf[i].size(); j++) {
                        smf_vec[j].push_back(vmf[i][j]);
                }
        }
    return smf_vec;
}

// return number of samples and number of variants in encoding file
void get_dimensions(string encodedTXT, int* dimensions, char delim){
	
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
				if (delim == ' '){
					size_t start = 0;
					size_t end = 0;
					const char delim = ' ';
					while ((start = line.find_first_not_of(delim, end)) != std::string::npos){
						end = line.find(delim, start);
						num_variants ++;
					}
				}else{
					num_variants = line.length();;
			
				}
			}
			num_samples ++;	
		}
		dimensions[0] = num_samples;
		dimensions[1] = num_variants;
	}
}
