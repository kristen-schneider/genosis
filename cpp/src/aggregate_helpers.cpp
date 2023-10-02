#include <algorithm>
#include <iostream>
#include <sstream>
#include <map>
#include <vector>
#include <string>
#include <fstream>

#include "aggregate_helpers.h"

using namespace std;

vector<string> get_query_samples_list(
        string query_samples_list) {
    vector<string> query_samples;
    string line;
    ifstream file(query_samples_list);
    while (getline(file, line)) {
        query_samples.push_back(line);
    }
    return query_samples;
}

map<int, vector<int>> get_chromosome_segments(
        vector<string> knn_results_file) {
    map<int, vector<int>> chromosome_segments;
    // iterate through all files in knn_results_file
    for (auto file_i : knn_results_file) {
        // file name format is chrmX.segmentYY.txt
        string chrom = file_i.substr(file_i.find("chrm") + 4, file_i.find("segment") - 5);
        int chromosome = stoi(chrom);
        string seg = file_i.substr(file_i.find("segment") + 7, file_i.find(".knn") - 7);
        int segment = stoi(seg);
        // add segment to list of segments for chromosome
        chromosome_segments[chromosome].push_back(segment);
    }
    return chromosome_segments;
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
        map<int, vector<int>> chromosome_segments,
        map<string, map<int, map<string, vector<float>>>> & query_chromosome_match_ID_segments
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
                // else split line on tab and add match ID to map
            else {
                string match_ID = line.substr(0, line.find("\t"));
                // if line is empty, skip
                if (match_ID.empty()) {
                    continue;
                }
                float score = stof(line.substr(line.find("\t") + 1));
                // if match_ID exists, add match ID to map
                if (query_chromosome_match_ID_segments[query_ID_line][chromosome].find(match_ID) !=
                    query_chromosome_match_ID_segments[query_ID_line][chromosome].end()) {
                    query_chromosome_match_ID_segments[query_ID_line][chromosome][match_ID][segment] = (1/score);
                }
                    // else add match ID to map
                else {
                    // initialize vector of zeros
                    sort(chromosome_segments[chromosome].begin(), chromosome_segments[chromosome].end(), greater<int>());
                    int max_segment = chromosome_segments[chromosome][0];
                    vector<float> zeros(max_segment, -1);
                    // add score to vector
                    // if (score == 0){ score = 10; }
                    zeros[segment] = score;
                    // add match ID and vector to map
                    query_chromosome_match_ID_segments[query_ID_line][chromosome][match_ID] = zeros;
                }
            }
        }
        file.close();
    }
        // if file is not open, print error
    else {
        cout << "Unable to open file " << filename << endl;
    }
}

void write_query_output(
        map<int, vector<int>> chromosome_segments,
        map<string, map<int, map<string, vector<float>>>> query_chromosome_match_ID_segments,
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
            for (int segment : chromosome_segments[chromosome_num]) {
                file << segment << ",";
            }
            file << endl;
            // write match IDs
            for (auto const& match_ID : chromosome.second) {
                file << match_ID.first << ",";
                // write scores
                for (float score : match_ID.second) {
                    file << score << ",";
                }
                file << endl;
            }
        }
    }
}

vector<string> read_ss_results_files(
        string ss_results_txt){
    vector<string> ss_results_files;
    ifstream file(ss_results_txt);
    string line;
    if (file.is_open()) {
        while (getline(file, line)) {
            // if line doesn't contain "segment", skip
            if (line.find("segment") == string::npos) {
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

map<string, vector<float>> score_samples(
        string query_hap_chrom){
    map<string, vector<float>> match_scores;
    string line;
    ifstream file(query_hap_chrom);
    string header;

    // format
    // MatchID,1,0,1,1,0,0,0,1,0...
    while (getline(file, line)) {
        vector<string> line_vector;
        string matchID;
        vector<float> segment_vector;

        // skip header
        if (header.empty()) {
            header = line;
            continue;
        }

        // split line on ',' and fill line_vector
        stringstream ss(line);
        while (ss.good()) {
            string substr;
            // add each substring to line_vector
            getline(ss, substr, ',');
            line_vector.push_back(substr);
        }
        // get matchID
        matchID = line_vector[0];
        // get segment vector
        for (int i = 1; i < line_vector.size(); i++) {
            // try stoi
            try {
                segment_vector.push_back(stof(line_vector[i]));
            } catch (invalid_argument) {
                // if stof fails, skip
                continue;
            }
        }
    
        // svs score = sum of scores reported by svs
        int svs_penalty = 10. // if score = -1 (doesn't appear in top k), penalty
        float matchID_svs_score = 0.;
        for (auto segment : segment_vector) {
            // if segment appears (score > 0), add svs score
            if (segment != -1){
                matchID_svs_score += segment;
            }else{
                matchID_svs_score += svs_penalty;
            }
        }    
    
        // popcount score = sum of appearances in top k
        float matchID_popcount_score = 0;
        for (auto segment : segment_vector) {
            // if segment appears (score > 0), add 1
            if (segment > 0){
                matchID_popcount_score += 1;
            }
        }

        // longest shared segment score = length of longest shared segment
        int matchID_lss = 0;
        int current_lss = 0;
        for (auto segment : segment_vector) {
            if (segment > 0) {
                current_lss++;
                matchID_lss = current_lss;
            }
            else {
                if (current_lss > matchID_lss) {
                    matchID_lss = current_lss;
                }
                current_lss = 0;
            }
        }
        matchID_lss = max(matchID_lss, current_lss);

        // cumulative shared segment = identify sum of all shared segments
        vector<int> shared_segments;
        int current_segment = 0;
        for (auto segment : segment_vector) {
            if (segment > 0) {
                current_segment++;
            }
            else {
                if (current_segment >= 3) {
                    shared_segments.push_back(current_segment);
                }
                current_segment = 0;
            }
        }
        if (current_segment > 0) {
            shared_segments.push_back(current_segment);
        }
        int matchID_sss = 0;
        for (auto segment : shared_segments) {
            matchID_sss += segment;
        }


        // IBD score
        // find all segments of length 3 or greater and sum their lengths
        // if there is a gap of 1 or less between segments, combine them
        vector<int> all_segments;
        vector<int> ibd_segments;
        vector<int> ibd_squared_segments;
        int ibd_curr_segment = 0;
        int ibd_prev_segment = 0;
        int gap_size = 0;

        for (auto segment : segment_vector) {
            if (segment > 0) {
                ibd_curr_segment++;
                gap_size = 0;
            }
            else{
                gap_size++;
                if (ibd_curr_segment >= 3) {
                    all_segments.push_back(ibd_curr_segment);
                }
                ibd_curr_segment = 0;

                if (gap_size == 1){
                    if (all_segments.size() > 1) {
                        ibd_curr_segment = all_segments.back();
                        all_segments.pop_back();
                        ibd_prev_segment = all_segments.back();
                        all_segments.pop_back();
                        ibd_curr_segment += ibd_prev_segment;
                        all_segments.push_back(ibd_curr_segment);
                        ibd_curr_segment = 0;
                    }
                }
                if (gap_size > 1){
                    if (all_segments.size() == 1) {
                        ibd_curr_segment = all_segments.back();
                        all_segments.pop_back();
                        ibd_segments.push_back(ibd_curr_segment);
                        ibd_squared_segments.push_back(ibd_curr_segment * ibd_curr_segment);
                        ibd_curr_segment = 0;
                    }
                }
            }

        }
        if (all_segments.size() > 0) {
            if (ibd_curr_segment >= 3) {
                if (gap_size <= 1){
                    ibd_prev_segment = all_segments.back();
                    all_segments.pop_back();
                    ibd_curr_segment += ibd_prev_segment;
                    all_segments.push_back(ibd_curr_segment);
                    ibd_curr_segment = 0;
                }
                else {
                    all_segments.push_back(ibd_curr_segment);
                }
            }
            else{
                ibd_curr_segment = all_segments.back();
                all_segments.pop_back();
                ibd_segments.push_back(ibd_curr_segment);
                ibd_squared_segments.push_back(ibd_curr_segment * ibd_curr_segment);
            }
        }
        else if (ibd_curr_segment >= 3){
            ibd_segments.push_back(ibd_curr_segment);
            ibd_squared_segments.push_back(ibd_curr_segment * ibd_curr_segment);
        }

        if (all_segments.size() == 1) {
            ibd_segments.push_back(all_segments.back());
            ibd_squared_segments.push_back(all_segments.back() * all_segments.back());
        }

        float matchID_ibd = 0;
        for (auto segment : ibd_segments) {
            matchID_ibd += segment;
        }

        float matchID_ibd_squared = 0;
        for (auto segment : ibd_squared_segments) {
            matchID_ibd_squared += segment;
        }

        // add match ID and scores to map
        match_scores[matchID] = {matchID_svs_score, matchID_popcount_score, matchID_lss, matchID_sss, matchID_ibd, matchID_ibd_squared};

    }
    file.close();
    return match_scores;
}

void write_all_chromosomes(string out_file,
                           string query_ID,
                           map<int, vector<pair<string, vector<float>>>> chromosome_match_ID_scores){
    ofstream file(out_file);
    // write header
    file << query_ID << endl;
    file << "Chromosome,MatchID,svs_score,popcount,longest_shared_segment,sum_shared_segment,ibd_score,ibd_score_squared" << endl;
    // for each chromosome
    for (auto const& chromosome : chromosome_match_ID_scores) {
        int chromosome_num = chromosome.first;
        // for each match ID
        for (auto const& match_ID_score : chromosome.second) {
            string match_ID = match_ID_score.first;
            float popcount_score = match_ID_score.second[0];
            float lss_score = match_ID_score.second[1];
            float sss_score = match_ID_score.second[2];
            float ibd_score = match_ID_score.second[3];
            float ibd_score_squared = match_ID_score.second[4];
            file << chromosome_num << ","
            << match_ID << ","
            << popcount_score << ","
            << lss_score << ","
            << sss_score << ","
            << ibd_score << ","
            << ibd_score_squared << endl;
        }
    }
}
