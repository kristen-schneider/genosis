#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <map>
#include <vector>

#include "aggregate_helpers.h"


using namespace std;

int main(int argc, char* argv[]) {
    
    string sim_search_results_dir = argv[1];
    string sim_search_results_filelist = argv[2];
    string query_results_dir = argv[3];

    // get list of chromosomes
    int num_chromosomes = 22;
    vector<int> chromosomes;
    for (int i = 1; i <= num_chromosomes; i++) {
        chromosomes.push_back(i);
    }

    // chromosome: list of segments for each chromosome
    map<int, vector<int> > chromosome_segments;

    // queryID: chromosome: matchID: [segment binary string]
    map<string, map<int, map<string, vector<float> > > > query_chromosome_match_ID_segments;

    // read faiss results file to get all files needed
    vector<string> sim_search_results_files;
    sim_search_results_files = read_ss_results_files(sim_search_results_filelist);
    chromosome_segments = get_chromosome_segments(sim_search_results_files);
    
    // iterate through all files in the sim search results directory
    for (auto file_i : sim_search_results_files) {
        string filename = sim_search_results_dir + file_i;
        string line;
        ifstream file(filename);
        string query_ID_line;

        // get chromosome and segment from file name
        // file name format: chrmX.segmentY.txt
        string chromosome_segment = filename.substr(filename.find_last_of("/\\") + 1);
        chromosome_segment = chromosome_segment.substr(4, chromosome_segment.find_last_of(".") - 4);
        int chromosome = stoi(chromosome_segment.substr(0, chromosome_segment.find("segment") - 1));
        int segment = stoi(chromosome_segment.substr(chromosome_segment.find("segment") + 7));
        // add segment to list of segments for chromosome
        // chromosome_segments[chromosome].push_back(segment);
        //cout << file_i << " " << chromosome << " " << segment << " " << endl;

        // read file and build datastructure
        read_QCMS(filename,
                    chromosome,
                    segment,
                    chromosome_segments,
                    query_chromosome_match_ID_segments);
    }
  
    // for all chromosomes in chromosome_segments, put segments in order
    map<int, vector<int> > sorted_chromosome_segments;
    for (auto chromosome : chromosome_segments) {
        // sort segments
        sort(chromosome.second.begin(), chromosome.second.end());
        sorted_chromosome_segments[chromosome.first] = chromosome.second;
    }
    
    // for each query, write out file
    write_query_output(sorted_chromosome_segments,
                        query_chromosome_match_ID_segments,
                        query_results_dir);

    // write output file segments.done
    string segments_done = query_results_dir + "segment_results.done";
    ofstream done_file(segments_done);
    done_file << query_chromosome_match_ID_segments.size();
}
