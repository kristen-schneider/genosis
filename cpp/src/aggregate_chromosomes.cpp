#include <algorithm>
#include <string>
#include <fstream>
#include <iostream>
#include <map>
#include <vector>

#include "aggregate_helpers.h"

using namespace std;

int main(int argc, char* argv[]) {

    string query_results_dir = argv[1];
    string query_samples_list = argv[2];

    // get list of chromosomes
    int num_chromosomes = 22;
    vector<int> chromosomes;
    for (int i = 1; i <= num_chromosomes; i++) {
        chromosomes.push_back(i);
    }
    
    // get list of queries
    vector<string> query_samples;
    query_samples = get_query_samples_list(query_samples_list);

    // full data structure
    // chromosome matchid: score
    map<int, vector<pair<string, vector<int>>>> chromosome_match_ID_scores;

    // for each query
    for (auto query : query_samples) {
        cout << "reading query: " << query << "\n";

        map<string, vector<int>> chrm_scores;
        string query_0 = query + "_0";
        string query_1 = query + "_1";

        // open directory for hap 0
        string query_hap0_dir = query_results_dir + query + "_0/";
        cout << "--hap: " << query_hap0_dir << "\n";
        // open file for each chromosome
        for (auto chromosome : chromosomes){
            string query_hap0_chrom = query_hap0_dir + "chrm" + to_string(chromosome) + ".csv";

            // if file exists, open it and read in all lines and store in map
            if (ifstream(query_hap0_chrom)) {
                chrm_scores = score_samples(query_hap0_chrom);
            }else {continue;}

            // for each match ID, add popcount and longest shared segment scores to map
            for (auto match : chrm_scores) {
            // get match_ID
            string match_ID = match.first;
            // get scores
            vector<int> scores = match.second;
            // add scores to map
            chromosome_match_ID_scores[chromosome].push_back(make_pair(match_ID,
                                                vector<int>{scores}));
            }
        }
        
        // all top matches have been collected for this query
        // write by chromosome to a file
        string out_file_0 = query_hap0_dir + query_0 + "_all_chromosomes.csv";
        write_all_chromosomes(out_file_0,
                              query_0,
                              chromosome_match_ID_scores);
        // clear chrm_scores
        chrm_scores.clear();
        // clear map
        chromosome_match_ID_scores.clear();

        // open directory for hap 1
        string query_hap1_dir = query_results_dir + query + "_1/";
        cout << "--hap: " << query_hap1_dir << "\n";
        // open file for each chromosome
        for (auto chromosome : chromosomes){
            
            string query_hap1_chrom = query_hap1_dir + "chrm" + to_string(chromosome) + ".csv";
            
            // if file exists, open it and read in all lines and store in map
            if (ifstream(query_hap1_chrom)) {
                chrm_scores = score_samples(query_hap1_chrom);
            }else {continue;}
            
            // for each match ID, add popcount and longest shared segment scores to map
            for (auto match : chrm_scores) {
            // get match_ID
            string match_ID = match.first;
            // get scores
            vector<int> scores = match.second;
            // add scores to map
            chromosome_match_ID_scores[chromosome].push_back(make_pair(match_ID,
                                                vector<int>{scores}));
            }
        }

        string out_file_1 = query_hap1_dir + query_1 + "_all_chromosomes.csv";
        write_all_chromosomes(out_file_1,
                              query_1,
                              chromosome_match_ID_scores);
        // clear chrm_scores
        chrm_scores.clear();
        // clear map
        chromosome_match_ID_scores.clear();
    }
}
