#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iterator>
#include <map>

#include "interpolate_map.h"

using namespace std;

int main(int argc, char *argv[]){
    string vcf_bps = argv[1];   // arg1 = vcf basepair file
    string map_file = argv[2];  // arg2 = plink map file
    string outfile = argv[3];   // arg3 = output file

    // read in vcf basepair file
    cout << "Reading in vcf basepair file: " << vcf_bps << endl;
    map<string, vector<int>> chr_bp_map = make_chr_bp_map(vcf_bps);

    // read in bp positions from plink map file and map to cm positions
    cout << "Reading in plink map file: " << map_file << endl;
    map<int, vector<tuple<int, float>>> chrm_cm_map = make_chr_cm_map(map_file);

    // interpolate map and write new output file
    cout << "Writing output file: " << outfile << endl;
    interpolate_map(chr_bp_map, chrm_cm_map, outfile);

    return 0;
}

map<string, vector<int>> make_chr_bp_map(string vcf_bps){
    // try to open file and exit if it fails
    ifstream file(vcf_bps);
    if (!file.is_open()){
        cout << "Error opening file: " << vcf_bps << endl;
        exit(1);
    }

    // read file line by line, split by whitespace, and store in vector
    string line;
    vector<int> basepair_positions;
    map<string, vector<int>> chr_bp_map;
    while (getline(file, line)){
        istringstream iss(line);
        vector<string> tokens{istream_iterator<string>{iss}, istream_iterator<string>{}};
        string chr = tokens[0];
        int bp = stoi(tokens[1]);
        chr_bp_map[chr].push_back(bp);
    }
    return chr_bp_map;
}

map<int, vector<tuple<int, float>>> make_chr_cm_map(string map_file){
    // try to open file and exit if it fails
    ifstream file(map_file);
    if (!file.is_open()){
        cout << "Error opening file: " << map_file << endl;
        exit(1);
    }

    // read file line by line, split by whitespace, and store in vector
    string line;
    int curr_chrom = 0;
    int chrm = 0;
    map<int, vector<tuple<int, float>>> chrm_cm_map;
    vector<tuple<int, float>> bp_cm_vectors;
    tuple<int, float> bp_cm;

    while (getline(file, line)){
        istringstream iss(line);
        vector<string> tokens{istream_iterator<string>{iss}, istream_iterator<string>{}};
        chrm = stoi(tokens[0]);
        float cm = stof(tokens[2]);
        int bp = stoi(tokens[3]);
        bp_cm = make_tuple(bp, cm);

        if (curr_chrom == 0){
            curr_chrom = chrm;
        }
        if (chrm != curr_chrom){
            // add bp_cm_map to chrm_cm_map
            chrm_cm_map[curr_chrom] = bp_cm_vectors;
            // reset bp_cm_vectors
            bp_cm_vectors.clear();
            // update curr_chrom
            curr_chrom = chrm;
            bp_cm_vectors.push_back(bp_cm);

        }else{
            bp_cm_vectors.push_back(bp_cm);
        }
    }
    // add bp_cm_map to chrm_cm_map at end of file
    chrm_cm_map[curr_chrom] = bp_cm_vectors;
    return chrm_cm_map;
}

void interpolate_map(
        map<string, vector<int>> chr_bp_map,
        map<int, vector<tuple<int, float>>> chrm_cm_map,
        string outfile) {

    // open output file to write to
    ofstream interpolated_map(outfile);
    if (!interpolated_map.is_open()) {
        cout << "Error opening file: " << outfile << endl;
        exit(1);
    }

    // iterate through each chromosome
    for (auto const &chr: chr_bp_map) {
        int vcf_bp_idx = 0;
        int map_bp_idx = 0;
        int curr_pos = -1;
        int map_file_bp = -1;
        // iterate through each basepair in the chromosome
        // convert chr.first ("chr8") to "8"
        int curr_chrm_int = stoi(chr.first.substr(3, chr.first.length()));for (auto const &bp: chr.second) {
            // while we are still in the current chromosome
            while (vcf_bp_idx < chr_bp_map[chr.first].size()) {
                curr_pos = chr_bp_map[chr.first][vcf_bp_idx];
                // get cm position from map file
                // get curr chrm map
                vector<tuple<int, float>> curr_chrm_map = chrm_cm_map[curr_chrm_int];
                // get current bp position from map file at map_bp_idx
                map_file_bp = get<0>(curr_chrm_map[map_bp_idx]);
                // if there is a bp entry in the map file
                if (curr_pos == map_file_bp) {
                    // the site has been mapped in reference, write to file
                    // chrm cm bp
                    interpolated_map << curr_chrm_int << " " << get<1>(curr_chrm_map[map_bp_idx]) << " "  << curr_pos << endl;
                    vcf_bp_idx++;
                    map_bp_idx++;
                // if the current bp is greater than the current map bp
                }
                // if the current bp is less than the current map bp
                else if (curr_pos < map_file_bp) {
                    // the site has not been mapped in reference, write to file
                    if (map_bp_idx == 0) {
                        // if the current bp is less than the first bp in the map file
                        // chrm cm bp
                        interpolated_map << curr_chrm_int << " " << get<1>(curr_chrm_map[map_bp_idx]) << " "  << curr_pos << endl;
                    } else {
                        // if the current bp is between two bp in the map file, interpolate
                        // get the previous bp and cm from the map file
                        int prev_map_bp = get<0>(curr_chrm_map[map_bp_idx - 1]);
                        float prev_map_cm = get<1>(curr_chrm_map[map_bp_idx - 1]);
                        // iterpolate
                        float frac = (float) (curr_pos - prev_map_bp) / (float) (map_file_bp - prev_map_bp);
                        float tmp_cm = prev_map_cm + frac * (get<1>(curr_chrm_map[map_bp_idx]) - prev_map_cm);
                        // chrm cm bp
                        interpolated_map << curr_chrm_int << " " << tmp_cm << " " << curr_pos << endl;
                    }
                    vcf_bp_idx++;
                } else if (curr_pos > map_file_bp) {
                    // if the current bp is greater than the current map bp
                    if (map_bp_idx == chrm_cm_map[curr_chrm_int].size() - 1) {
                        // if the current bp is greater than the last bp in the map file
                        // chrm cm bp
                        interpolated_map << curr_chrm_int << " " << get<1>(curr_chrm_map[map_bp_idx]) << " "  << curr_pos << endl;
                        vcf_bp_idx++;
                    } else {
                        // increment the map_bp_idx
                        map_bp_idx++;
                    }
                }
            }
        }
    }
}
