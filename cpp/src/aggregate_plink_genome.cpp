#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <map>
#include <vector>

#include "aggregate_plink_genome.h"

using namespace std;

int main(int argc, char* argv[]){
    string plink_file = argv[1];
    int k = stoi(argv[2]);
    string output_file = argv[3];
    string query_file = argv[4];

    vector<string> queries;
    char delim = ' ';

    // create plink map from plink genome file
    map<string, map<string, float>> plink_map;
    plink_map = main_plink_aggregate(plink_file, delim);

    // sort plink map by value for each sample
    map<string, vector<pair<string, float>>> sorted_plink_map;
    sorted_plink_map = sort_full_map(plink_map);

    // return top k for each query sample
    queries = read_query_file(query_file);
    map<string, vector<pair<string, float>>> top_k_map;
    top_k_map = return_top_k(sorted_plink_map, queries, k);

    // write top k to file
    write_top_k(top_k_map, k, output_file);

    cout << "DONE" << endl;
}

map<string, map<string, float>> main_plink_aggregate(string plink_file, char delim=' '){
    /*
     * open sample-wise plink file
     */
    // TO RETURN:
    // map to store the sampleID_1 sampleID_2 dist
    map<string, map<string, float>> plink_map;
    ifstream plink_file_stream;
    plink_file_stream.open(plink_file);
    cout << "READING PLINK GENOME FILE..." << plink_file << endl;
    if (!plink_file_stream.is_open()){
        cout << "FAILED TO OPEN: " << plink_file << endl;
        exit(1);
    }
    string header = "";
    string line;
    while (getline(plink_file_stream, line)) {
        if (header == "") {
            header = line;
        } else {
            // map to store sampleID and dist
            map<string, float> sampleID1_dist;
            map<string, float> sampleID2_dist;

            // split line and assign value by column
            vector<string> vec_line;
            split_line(line, delim, vec_line);
            string sampleID_1 = vec_line[1];
            string sampleID_2 = vec_line[3];
            float dist_value = stof(vec_line[11]);

            // fill map with sample 1 as key
            if (plink_map.find(sampleID_1) != plink_map.end()) {
                plink_map[sampleID_1][sampleID_2] = dist_value;
            } else {
                sampleID1_dist[sampleID_2] = dist_value;
                plink_map[sampleID_1] = sampleID1_dist;
            }
            // fill map with sample 2 as key
            if (plink_map.find(sampleID_2) != plink_map.end()) {
                plink_map[sampleID_2][sampleID_1] = dist_value;
            } else {
                sampleID2_dist[sampleID_1] = dist_value;
                plink_map[sampleID_2] = sampleID2_dist;
            }
        }
    }
    cout << "DONE READING PLINK GENOME FILE" << endl;
    return plink_map;
}

vector<pair<string, float>> sort_map(map<string, float> str_flt_map){
    /*
     * returns a string-float map sorted by float value from largest to smallest
     */
    vector<pair<string, float>> vec;
    copy(str_flt_map.begin(), str_flt_map.end(), back_inserter<vector<pair<string, float>>>(vec));
    sort(vec.begin(), vec.end(), [](const pair<string, float> &l, const pair<string, float> &r) {
        return l.second > r.second;
    });
    return vec;
}

// sort full map for each sample
map<string, vector<pair<string, float>>> sort_full_map(map<string, map<string, float>> str_str_flt_map){
    /*
     * returns a map of maps sorted by float value from largest to smallest
     */
    cout << "SORTING FULL MAP" << endl;
    map<string, vector<pair<string, float>>> sorted_map;
    for (auto &p : str_str_flt_map) {
//        cout << "SORTING: " << p.first << endl;
        vector<pair<string, float>> sorted_vec;
        sorted_vec = sort_map(p.second);
        sorted_map[p.first] = sorted_vec;
    }
    cout << "DONE SORTING FULL MAP" << endl;
    return sorted_map;
}


// return top k for each sample
map<string, vector<pair<string, float>>> return_top_k(map<string, vector<pair<string, float>>> str_vec_map, vector<string> queries, int k){
    /*
     * returns the top k values for each sample in a vector of queries
     */
    cout << "COMPUTING TOP K" << endl;
    map<string, vector<pair<string, float>>> top_k_map;
    for (auto query : queries) {
        vector<pair<string, float>> top_k;
        int i = 0;
        for (auto &p : str_vec_map[query]) {
            if (i < k) {
                top_k.push_back(p);
                i++;
            } else {
                break;
            }
        }
        top_k_map[query] = top_k;
    }
    return top_k_map;
}


void write_top_k(map<string, vector<pair<string, float>>> top_k_map, int k, string output_file){
    /*
     * writes the top k values for each sample to a file
     */
    cout << "WRITING TOP K" << endl;
    ofstream output_file_stream;
    output_file_stream.open(output_file);
    if (!output_file_stream.is_open()){
        cout << "FAILED TO OPEN: " << output_file << endl;
        exit(1);
    }
    for (auto &p : top_k_map) {
        output_file_stream << p.first << "\t";
        for (auto &q : p.second) {
            output_file_stream << q.first << "\t" << q.second << "\t";
        }
        output_file_stream << endl;
    }
    output_file_stream.close();
    cout << "DONE WRITING TOP K" << endl;
}

void split_line(const string &s, char delim, vector<string> &elems){
    /*
     * Splits a string by white space
     */
    stringstream ss;
    ss.str(s);
    string item;
    while (getline(ss, item, delim)) {
        if (item != "") {
            elems.push_back(item);
        }
    }
}

// read file of queries
vector<string> read_query_file(string query_file){
    /*
     * reads a file of queries and returns a vector of queries
     */
    vector<string> query_vec;
    ifstream query_file_stream;
    query_file_stream.open(query_file);
    if (!query_file_stream.is_open()){
        cout << "FAILED TO OPEN: " << query_file << endl;
        exit(1);
    }
    string line;
    while (getline(query_file_stream, line)) {
        query_vec.push_back(line);
    }
    return query_vec;
}
