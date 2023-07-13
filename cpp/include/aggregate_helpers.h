#ifndef AGGREGATE_HELPERS_H
#define AGGREGATE_HELPERS_H

#endif //AGGREGATE_HELPERS_H


#include <string>
#include <map>
#include <vector>

using namespace std;

vector<pair<string, int>> read_match_file(
        string query_hap_chrom);

vector<pair<string, int>> sort_matches(
        vector<pair<string, int>> match_scores);

vector<string> get_query_samples_list(
        string query_samples_list);

vector<string> read_ss_results_files(
        string ss_results_txt);

void read_QCMS(
        string filename,
        int chromosome,
        int segment,
        map<string, map<int, map<string, vector<int>>>> & query_chromosome_match_ID_segments
);

void write_query_output(
        map<int, vector<int>> chromosome_segments,
        map<string, map<int, map<string, vector<int>>>> query_chromosome_match_ID_segments,
        string query_results_dir
);

void write_all_chromosomes(string out_file,
                          string query_ID,
                          map<int, vector<pair<string, int>>> chromosome_match_ID_scores);
