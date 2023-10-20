//
// Created by Kristen Schneider on 7/5/23.
//

#ifndef AGGREGATE_AGGREGATE_H
#define AGGREGATE_AGGREGATE_H

#endif //AGGREGATE_AGGREGATE_H

#include <string>
#include <map>
#include <vector>

using namespace std;

map<int, vector<int> > get_chromosome_segments(
        vector<string> chromosome_segments_file);

map<string, vector<float> > score_samples(
        string query_hap_chrom);

vector<string> get_query_samples_list(
        string query_samples_list);

vector<string> read_ss_results_files(
        string ss_results_txt);

void read_QCMS(
        string filename,
        int chromosome,
        int segment,
        map<int, vector<int> > chromosome_segments,
        map<string, map<int, map<string, vector<float> > > > & query_chromosome_match_ID_segments
);

void write_query_output(
        map<int, vector<int> > chromosome_segments,
        map<string, map<int, map<string, vector<float> > > > query_chromosome_match_ID_segments,
        string query_results_dir
);

void write_all_chromosomes(string out_file,
                          string query_ID,
                          map<int, vector<pair<string, vector<float> > > > chromosome_match_ID_scores);
