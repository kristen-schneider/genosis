#ifndef ENCODEPOS_H
#define ENCODEPOS_H

#endif //ENCODEPOS_H

#include <cstdlib>
#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <htslib/vcfutils.h>
#include <htslib/kstring.h>
#include <htslib/kseq.h>
#include <htslib/synced_bcf_reader.h>
#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <vector>

using namespace std;

void encode_positions(string map_file,
		string vcf_slice_file,
		string sample_IDs_file, 
		string output_encoding_file);
void write_SMF(vector<string> all_sample_IDs,
                vector<vector<float>> smf,
                string output_encoding_file);
