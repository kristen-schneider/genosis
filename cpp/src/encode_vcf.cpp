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
	cout << "Encoding vectors..." << endl;
	// read config file
	cout << "...Reading config File..." << endl;
	string configFile = argv[1];   // configuration file will all options
	map<string, string> config_options;
	config_options = get_config_options(configFile);

	string sample_IDs_file = argv[2];
	string vcf_slice_file = argv[3];
	string output_encoding_file = argv[4];

	// access each option by variable name
	string encoding_file = config_options["encoding_file"];
	int slice_size = stoi(config_options["slice_size"]);
	string out_dir = config_options["out_dir"];
    	string output_base_name = config_options["out_base_name"];
	cout << "...Done reading config file" << endl;

	// make encoding map
	//map<string, int> encoding_map = make_encoding_map(encoding_file);
	cout << "...Making encoding map..." << endl;
	map<string, vector<int>> biallelic_encoding_map = make_biallelic_encoding_map(encoding_file); 
	cout << "...Done making encoding map." << endl;

	//encode_vcf(sample_IDs_file, vcf_slice_file, encoding_map, output_encoding_file);
	cout << "...Ecoding VCF..." << endl;
	encode_vcf(sample_IDs_file, vcf_slice_file, biallelic_encoding_map, output_encoding_file);
	cout << "Done." << endl;
	return 0;
}

/*
 * takes one sliced vcf file and encoding instructions
 * and writes output encoiding file
 */
void encode_vcf(string sample_IDs_file, string vcf_slice_file, map<string, vector<int>> encoding_map, string output_encoding_file){
	
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
	
	// initialize and allocate bcf1_t object
        bcf1_t *vcf_record = bcf_init();
        if (vcf_record == NULL) {
        	cout << "ERROR: record is empty" << endl;
        }

	//vector<vector<int>> all_genotype_encodings; // vector of vectors
	vector<vector<int>> all_haplotype_encodings;
	//vector<vector<int>> all_haplotypeA_encodings; // vector of vectors
	//vector<vector<int>> all_haplotypeB_encodings; // vector of vectors


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
		//vector<int> record_encoding_vector;
		vector<int> haplotype_encoding_vector;
		//vector<int> haplotypeA_encoding_vector;
		//vector<int> haplotypeB_encoding_vector;

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
			//record_encoding_vector.push_back(encoding_map[s_gt]);
			vector<int> biallelic_encoding = encoding_map[s_gt];
			haplotype_encoding_vector.push_back(biallelic_encoding[0]);
			haplotype_encoding_vector.push_back(biallelic_encoding[1]);
			//haplotypeA_encoding_vector.push_back(biallelic_encoding[0])
			//haplotypeB_encoding_vector.push_back(biallelic_encoding[1])
				
		}
		//all_genotype_encodings.push_back(record_encoding_vector);
		all_haplotype_encodings.push_back(haplotype_encoding_vector);
		//all_haplotypeA_encodings.push_back(haplotypeA_encoding_vector);
		//all_haplotypeB_encodings.push_back(haplotypeB_encoding_vector);

		//record_encoding_vector.clear();
		haplotype_encoding_vector.clear();
		//haplotypeA_encoding_vector.clear();
		//haplotypeB_encoding_vector.clear();
		
	} // end of reading records
	// transposing data
	cout << "...transposing data..." << endl;
	//vector<vector<int>> sample_major_format_vec = transpose(all_genotype_encodings);
	vector<vector<int>> sample_major_format_hap_vec = transpose(all_haplotype_encodings); 
	//vector<vector<int>> sample_major_format_hapA = transpose(all_haplotypeA_encodings);
	//vector<vector<int>> sample_major_format_hapB = transpose(all_haplotypeB_encodings);

	vector<string> all_sample_IDs = get_sample_IDs(sample_IDs_file); 
	// writing smf
	cout << "...writing sample major format encodings to file..." << endl;
	//write_SMF(all_sample_IDs, sample_major_format_vec, output_encoding_file);	
	cout << sample_major_format_hap_vec.size() << endl;
	write_SMF(all_sample_IDs, sample_major_format_hap_vec, output_encoding_file);
	//write_SMF_haplotype(all_sample_IDs, sample_major_format_hapA, sample_major_format_hapB, output_encoding_file);	
}	

/*
 * given a set of sample IDs,
 * and a set of smf encoded vectors, 
 * write: sampleID encoding
 */
void write_SMF(vector<string> all_sample_IDs, vector<vector<int>> smf, string output_encoding_file){
	
	// open output file to write encoding
	ofstream output_stream;
	output_stream.open(output_encoding_file);

	// format sample ID float float float...
	// space delim
	int SID_i = 0;
	int binary = -1;
	for (int i = 0; i < smf.size(); i++) {
		if (binary == -1){
			binary = 0;
		}else if (binary == 0){
			binary = 1;		
		}else{
			SID_i ++;
			binary = 0;
		}
		vector<int> sample = smf.at(i);
		output_stream << all_sample_IDs.at(SID_i) << " ";
		for(int j = 0; j < sample.size(); j++) {
			output_stream << sample.at(j) << " ";
		}
		output_stream << endl;
	}
}
