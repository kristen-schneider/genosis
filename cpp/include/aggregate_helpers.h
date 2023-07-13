#ifndef AGGREGATE_HELPERS_H
#define AGGREGATE_HELPERS_H

#endif //AGGREGATE_HELPERS_H

#include <string>
#include <map>
#include <vector>

using namespace std;

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
