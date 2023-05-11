#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iterator>
#include <map>
#include <tuple>

# include "segment_boundary_map.h"

using namespace std;


int main(int argc, char *argv[]){
    string interpolated_map = argv[1];
    string outfile = argv[2];

    float slice_size = 1;

    map<int, vector<tuple<int, int>>> slice_boundary_map = make_slice_boundary_map(interpolated_map, slice_size);
    write_slice_boundary_map(slice_boundary_map, outfile);

    return 0;
}

map<int, vector<tuple<int, int>>> make_slice_boundary_map(
        string interpolated_map,
        float slice_size){
    // try to open file and exit if it fails
    ifstream file(interpolated_map);
    if (!file.is_open()){
        cout << "Error opening file: " << interpolated_map << endl;
        exit(1);
    }

    // read file line by line, split by whitespace, and store in vector
    string line;
    map<int, vector<tuple<int, int>>> slice_boundary_map;
    vector<tuple<int, int>> chrm_start_end;

    float curr_slize_size = 0;
    int start = 0;
    int end = 0;
    int curr_chrm = 0;
    int segment_idx = 0;

    while (getline(file, line)){
        istringstream iss(line);
        vector<string> tokens{istream_iterator<string>{iss}, istream_iterator<string>{}};

        int chrm = stoi(tokens[0]);
        float cm = stof(tokens[1]);
        int bp = stoi(tokens[2]);

        if (curr_chrm == 0){
            curr_chrm = chrm;
        }
        // onto the next chromosome
        if (chrm != curr_chrm){
            // add last segment from previous chromosome
            chrm_start_end.emplace_back(start, end);
            slice_boundary_map[curr_chrm] = chrm_start_end;
            curr_chrm = chrm;
            chrm_start_end.clear();
            segment_idx = 0;
            curr_slize_size = 0;
            start = 0;
            end = 0;
        }
        // on the same chromosome
        else{
            if (curr_slize_size < slice_size + segment_idx){
                curr_slize_size = cm;
		end = bp;
            }
            else{
                end = bp;
                chrm_start_end.emplace_back(start, end);
                start = end;
                segment_idx++;
                curr_slize_size = segment_idx;
            }
        }
    }
    // add the last segment
    chrm_start_end.emplace_back(start, end);
    slice_boundary_map[curr_chrm] = chrm_start_end;
    return slice_boundary_map;
}

void write_slice_boundary_map(
        map<int, vector<tuple<int, int>>> slice_boundary_map,
        string outfile){
    ofstream file(outfile);
    if (!file.is_open()){
        cout << "Error opening file: " << outfile << endl;
        exit(1);
    }
    int segment_idx = 0;
    for (auto const& [chrm, start_end] : slice_boundary_map){
        segment_idx = 0;       
        for (auto const& [start, end] : start_end){
            file << chrm << "\t" << segment_idx << "\t" << start << "\t" << end << endl;
            segment_idx++;
        }
    }
}

