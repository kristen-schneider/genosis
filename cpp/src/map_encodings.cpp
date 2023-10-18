#include <cstdlib>
#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <vector>
#include <sstream>
#include <iterator>

#include "map_encodings.h"

using namespace std;

map<string, vector<int> > map_gt_encoding(string encoding_file){
    ifstream file(encoding_file);
    if (!file.is_open()){
        cout << "Error opening file: " << encoding_file << endl;
        exit(1);
    }

    string line;
    map<string, vector<int> > gt_encoding_map;
    vector<int> single_gt_encoding;

    while (getline(file, line)){
        istringstream iss(line);
        vector<string> tokens{istream_iterator<string>{iss}, istream_iterator<string>{}};

        string gt = tokens[0];
        int encoding_1 = stoi(tokens[1]);
        int encoding_2 = stoi(tokens[2]);

        single_gt_encoding.push_back(encoding_1);
        single_gt_encoding.push_back(encoding_2);

        gt_encoding_map[gt] = single_gt_encoding;
        single_gt_encoding.clear();
    }
    return gt_encoding_map;
}
