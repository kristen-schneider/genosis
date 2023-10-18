#include <cstdlib>
#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <vector>
#include <sstream>
#include <iterator>

#include "read_map.h"

using namespace std;

/**
 * creates a map from basepair position to centimorgan position for all chromosomes
 *
 * @param interpolated_map_file MAP file with format: chrm cm bp
 *
 * @return chrm_bp_cm_map chrm: bp: cm
 */
map<int, map<int, float> > map_bp_cm(string interpolated_map_file){
    ifstream file(interpolated_map_file);
    if (!file.is_open()){
        cout << "Error opening file: " << interpolated_map_file << endl;
        exit(1);
    }

    int curr_chrm = 0;

    string line;
    map<int, map<int, float> > chrm_bp_cm_map;
    map<int, float> single_bp_cm;

    while (getline(file, line)){
        istringstream iss(line);
        vector<string> tokens{istream_iterator<string>{iss}, istream_iterator<string>{}};

        int chrm = stoi(tokens[0]);
        int bp = stoi(tokens[2]);
        float cm = stof(tokens[1]);

        if (curr_chrm == 0){
            curr_chrm = chrm;
            single_bp_cm[bp] = cm;
        }
        else if (chrm != curr_chrm){
            chrm_bp_cm_map[curr_chrm] = single_bp_cm;
            curr_chrm = chrm;
            single_bp_cm.clear();
            single_bp_cm[bp] = cm;
        }
        else{
            single_bp_cm[bp] = cm;
        }

    }
    chrm_bp_cm_map[curr_chrm] = single_bp_cm;
    return chrm_bp_cm_map;
}
