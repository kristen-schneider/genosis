#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <map>
#include <vector>

#include "write_query_results.h"

using namespace std;

/*
 * Read a file which lists all files in a directory
 * @param faiss_results_txt: file will all faiss results in it
 * @return faiss_results_files: vector of all files with faiss results
 */
vector<string> read_ss_results_files(
        string ss_results_txt){
    vector<string> ss_results_files;
    ifstream file(ss_results_txt);
    string line;
    if (file.is_open()) {
        while (getline(file, line)) {
            // if this is not a results file, dont add it to the list
            if (line.find("knn") == string::npos){
                continue;
            }
            ss_results_files.push_back(line);
        }
        file.close();
    }
    else {
        cout << "Unable to open file " << ss_results_txt << endl;
    }
    return ss_results_files;
}

/*
 * Read knn file and add match IDs to map
 * @param filename: name of file to read
 * @param chromosome: chromosome number
 * @param segment: segment number
 * @param query_chromosome_match_ID_segments: map of query ID to chromosome to match IDs to segments
 * @return: void
*/
void read_QCMS(
        string filename,
        int chromosome,
        int segment,
        map<string, map<int, map<string, vector<int>>>> & query_chromosome_match_ID_segments
){
    ifstream file(filename);
    string line;
    string query_ID_line;

    // if file is open, read file line by line
    if (file.is_open()) {
    
        // read file line by line.
        // format:
        // Query: queryID
        // matchID distance
        // matchID distance
        //
        // Query: queryID
        // matchID distance
        // matchID distance
        // ...
    
        while (getline(file, line)) {
        // if "Query: " is found, get query ID
            if (line.find("Query: ") != string::npos) {
                query_ID_line = line;
                query_ID_line = query_ID_line.substr(7);
            }
            // else split the line on tab and add match ID to map
            else {
                string match_ID = line.substr(0, line.find("\t"));
                // if line is empty (last one), skip
                if (match_ID.empty()) {
                    continue;
                }
            
                // if match ID exists, add segment to match ID in map
                if (query_chromosome_match_ID_segments[query_ID_line][chromosome].find(match_ID) !=
                    query_chromosome_match_ID_segments[query_ID_line][chromosome].end()) {
                    query_chromosome_match_ID_segments[query_ID_line][chromosome][match_ID].push_back(segment);
                }
            
                // else add match ID to map
                else{
                    query_chromosome_match_ID_segments[query_ID_line][chromosome][match_ID] = {segment};
                }         
            }
        }
        file.close();
    }
    // if file did not open, print error
    else{
        cout << "Unable to open file " << filename << endl;
    }
}


void write_query_output(
        map<int, vector<int>> chromosome_segments,
        map<string, map<int, map<string, vector<int>>>> query_chromosome_match_ID_segments,
        string query_results_dir
        ){
    // for each query make a directory for output
    for (auto const& query : query_chromosome_match_ID_segments) {
        string query_ID = query.first;
        // make directory for this query if it doesn't exist
        
        string query_dir = query_results_dir + query_ID;
        string mkdir_command = "mkdir " + query_dir;
        system(mkdir_command.c_str());

        // for each chromosome write out file
        // format:
        // segment, 0, 1, 2, 3, ...
        // matchID1, 1, 0, 0, 1, ...
        // matchID2, 0, 1, 1, 0, ...
        // ...

        for (auto const& chromosome : query.second) {
            int chromosome_num = chromosome.first;
            string query_chrm_file = query_dir + "/chrm" + to_string(chromosome_num) + ".csv";
            ofstream file(query_chrm_file);
            // write header
            file << "segment,";
            int num_segments = chromosome_segments[chromosome_num].size();
            for (int i = 0; i < num_segments; i++) {
                file << i << ",";
            }
            file << endl;
            // write match IDs
            for (auto const& match_ID : chromosome.second) {
                file << match_ID.first << ",";
                // for each segment in chromosome,
                // write 1 if match ID is in segment,
                // 0 otherwise
                for (int segment : chromosome_segments[chromosome_num]) {
                    if (find(match_ID.second.begin(), match_ID.second.end(), segment) != match_ID.second.end()) {
                        file << "1,";
                    }
                    else {
                        file << "0,";
                    }
                }
                file << endl;
            }
        }
    }
}

