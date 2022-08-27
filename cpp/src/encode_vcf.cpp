#include <iostream>
#include <fstream>
#include <math.h>
#include <map>
#include <string>
#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <htslib/vcfutils.h>
#include <htslib/kstring.h>
#include <htslib/kseq.h>
#include <htslib/synced_bcf_reader.h>
#include <vector>

#include "encode_vcf.h"
#include "map_encodings.h"
#include "read_config.h"
#include "utils.h"

using namespace std;

int main(int argc, char* argv[]){
	// read config file
	cout << "Reading Config File..." << endl;
	string configFile = argv[1];   // configuration file will all options
	map<string, string> config_options;
	config_options = get_config_options(configFile);

	string vcf_slice_file = argv[2];
	string output_encoding_file = argv[3];

	// access each option by variable name
	string encoding_file = config_options["encoding_file"];
	int slice_size = stoi(config_options["slice_size"]);
	string out_dir = config_options["out_dir"];
    	string output_base_name = config_options["out_base_name"];

	// make encoding map
	map<string, int> encoding_map = make_encoding_map(encoding_file);
	
	encode_vcf(vcf_slice_file, encoding_map, output_encoding_file);

	return 0;
}

/*
 * takes one sliced vcf file and encoding instructions
 * and writes output encoiding file
 */
void encode_vcf(string vcf_slice_file, map<string, int> encoding_map, string output_encoding_file){
	
	// converts vcfFile name to const char for htslib
	const char *vcf_slice = vcf_slice_file.c_str();
	// open VCF file with htslib
	cout << "VCF; " << vcf_slice << endl;
	htsFile *vcf_stream = bcf_open(vcf_slice, "r");
	if (!vcf_stream){
		cout << "FAILED TO OPEN: " << vcf_slice << endl;
		exit(1);
	}
		
	bcf_hdr_t *vcf_header = bcf_hdr_read(vcf_stream);
	if (vcf_header == NULL) {
        	throw runtime_error("Unable to read header.");
		exit(1);
        }
	
	// getting number of samples
	int num_samples = get_num_samples(vcf_header);
	cout << "NUM SAMPLES: " << num_samples << endl;
	// gettting all sequence names
	const char **sequence_names = get_sequence_names(vcf_header);
	
	// initialize and allocate bcf1_t object
        bcf1_t *vcf_record = bcf_init();
        if (vcf_record == NULL) {
        	cout << "ERROR: record is empty" << endl;
        }

	vector<vector<int>> all_genotype_encodings; // vector of vectors

	cout << "...reading genotypes." << endl;
	while (bcf_read(vcf_stream, vcf_header, vcf_record) == 0){

		bcf_unpack(vcf_record, BCF_UN_ALL);
		bcf_unpack(vcf_record, BCF_UN_INFO);

		// making variable names for each column
		string chrm = bcf_hdr_id2name(vcf_header, vcf_record->rid);
		int pos = (unsigned long)vcf_record->pos;
		string id = vcf_record->d.id;
		string ref = vcf_record->d.allele[0];
		string alt = vcf_record->d.allele[1];
		double_t qual = vcf_record->qual;

		// make chromosome through quality into one vector of strings
		vector<string> chrm_thru_qual_vector = 
			{chrm, to_string(pos), id, ref, alt, to_string(qual)};
		string s;
		for (const auto &col : chrm_thru_qual_vector) s += "\t" + col;
		string chrm_thru_qual_string;
		chrm_thru_qual_string = chrm_thru_qual_string +
			chrm + "\t" +
			to_string(pos) + "\t" +
			id + "\t" +
			ref + "\t" +
			alt + "\t" +
			to_string(qual) + "\t";

		// reading genotypes
		string s_gt;
		int *gt = NULL;
		int num_alleles = 0;
		int ngt_arr = 0;
		vector<int> record_encoding_vector;

		num_alleles = bcf_get_genotypes(vcf_header, vcf_record,  &gt, &ngt_arr);
		int alleles_per_gt = num_alleles/num_samples;
		for (int i = 0; i < num_samples; i++){
			int allele1 = bcf_gt_allele(gt[i*alleles_per_gt+0]);
			int allele2;
			
			// replace unknowns
			if (allele1 == -1){ s_gt = ".|.";}
			else{
				allele2 = bcf_gt_allele(gt[i*alleles_per_gt+1]);
				s_gt = to_string(allele1)+"|"+to_string(allele2);
			}
			record_encoding_vector.push_back(encoding_map[s_gt]);
		}
		all_genotype_encodings.push_back(record_encoding_vector);
		record_encoding_vector.clear();
	} // end of reading records
	// transposing data
	cout << "...transposing data..." << endl;
	vector<vector<int>> sample_major_format_vec = transpose(all_genotype_encodings);
	// writing smf
	cout << "...writing sample major format encodings to file..." << endl;
	write_SMF(sample_major_format_vec, output_encoding_file);
	
}


int get_num_samples(bcf_hdr_t *vcf_header){
	/*
	 * returns number of samples in VCF file
	 */

	int num_samples = -1;
	num_samples = bcf_hdr_nsamples(vcf_header);
	return num_samples;
}

const char **get_sequence_names(bcf_hdr_t *vcf_header){
	/*
	 * returns sequnce IDs/names in VCF file
	 */
	int n = 0;
	
	const char **sequence_names = NULL;
	bcf_hdr_seqnames(vcf_header, &n);
	
	return sequence_names;
}

void write_SMF(vector<vector<int>> smf, string output_encoding_file){
	/*
	 * writes sample major format
	 * to outfile
	 */
	
	// open output file to write encoding
	ofstream output_stream;
	output_stream.open(output_encoding_file);

	for (int i = 0; i < smf.size(); i++) {
		vector<int> sample = smf.at(i);
		for(int j = 0; j < sample.size(); j++) {
			output_stream << sample.at(j);
		}
		output_stream << endl;
	}
}
