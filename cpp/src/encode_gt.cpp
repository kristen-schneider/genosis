#include <cstdlib>
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

#include "encode_gt.h"
#include "map_encodings.h"
#include "read_config.h"
#include "utils.h"

using namespace std;

/*
 * takes one sliced vcf file and encoding instructions
 * and writes output encoiding file
 */
void encode_gt_vectors(string sample_IDs_file, 
		string vcf_slice_file, 
		map<string, vector<int>> encoding_map, 
		string output_encoding_file){
	
	// converts vcfFile name to const char for htslib
	const char *vcf_slice = vcf_slice_file.c_str();
	
	// open VCF file with htslib
	htsFile *vcf_stream = bcf_open(vcf_slice, "r");
	if (!vcf_stream){
		cout << "FAILED TO OPEN: " << vcf_slice << endl;
		exit(1);
	}	

	// read vcf header
	bcf_hdr_t *vcf_header = bcf_hdr_read(vcf_stream);
	if (vcf_header == NULL) {
		throw runtime_error("Unable to read header.");
		exit(1);
        } 

	// get number of samples
	int num_samples;
	num_samples = bcf_hdr_nsamples(vcf_header);
	char ** sample_ids = vcf_header->samples;

	// allocate space for vcf record
        bcf1_t *vcf_record = bcf_init();
        if (vcf_record == NULL) {
        	cout << "ERROR: record is empty" << endl;
        }

	vector<vector<int>> all_haplotype_encodings;
	vector<int> all_bp_positions;
	
	cout << "...Reading genotypes." << endl;
	while (bcf_read(vcf_stream, vcf_header, vcf_record) == 0){
		bcf_unpack(vcf_record, BCF_UN_ALL);
		bcf_unpack(vcf_record, BCF_UN_INFO);

		// storing data for each column
		string chrm = bcf_hdr_id2name(vcf_header, vcf_record->rid);	// chrom
		int pos = (unsigned long)vcf_record->pos;			// pos
		string id = vcf_record->d.id;					// sample id
		string ref = vcf_record->d.allele[0];				// ref allele
		string alt = vcf_record->d.allele[1];				// alt allele
		double_t qual = vcf_record->qual;				// quality
		all_bp_positions.push_back(pos);				// add pos to vector
	
		// reading genotypes
		string s_gt;
		int *gt = NULL;
		int number_total_alleles = 0;
		int ngt_arr = 0;
		vector<int> haplotype_encoding_vector;
		number_total_alleles = bcf_get_genotypes(vcf_header, vcf_record,  &gt, &ngt_arr);
		int alleles_per_gt = number_total_alleles/num_samples;
		
		// for each sample in the record, convert gt to encoding
		for (int i = 0; i < num_samples; i++){

			// separate two alleles for each gt
			int allele1 = bcf_gt_allele(gt[i*alleles_per_gt+0]);
			int allele2 = bcf_gt_allele(gt[i*alleles_per_gt+1]);
			
			// replace unknowns and concatinate genotypes to " | " format
			if (allele1 != 0 and allele1 != 1 and allele1 != 2 and allele1 != 3){
				s_gt = ".|" + to_string(allele2);
			}
			if (allele2 != 0 and allele2 != 1 and allele2 != 2 and allele2 != 3){
				s_gt = to_string(allele1) + "|.";
                        }
			else{
				s_gt = to_string(allele1)+"|"+to_string(allele2);
			}

			// mapping genotype to encoding vector
			vector<int> biallelic_encoding = encoding_map[s_gt];
			haplotype_encoding_vector.push_back(biallelic_encoding[0]);
			haplotype_encoding_vector.push_back(biallelic_encoding[1]);
		}
		all_haplotype_encodings.push_back(haplotype_encoding_vector);
		haplotype_encoding_vector.clear();
	} // end of reading records
	cout << "...Done reading genotypes." << endl;
	
	// transposing data
	cout << "...Transposing data." << endl;
	vector<vector<int>> sample_major_format_hap_vec = transpose_int(all_haplotype_encodings); 
	cout << "...Done transposing data." << endl;

	vector<string> all_sample_IDs = get_sample_IDs(sample_IDs_file); 
	// writing smf
	cout << "...Writing sample major format encodings to file: " << output_encoding_file << "." << endl;
	write_SMF(all_sample_IDs, sample_major_format_hap_vec, output_encoding_file);
	cout << "...Done writing sample major format." << endl;

	// writing positinal encoding
	//cout << "...writing positional encodings to file..." << endl;
	//write_positional_encoding(all_bp_positions, all_sample_IDs, sample_major_format_hap_vec, output_position_file);	
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
		output_stream << all_sample_IDs[SID_i] << "_" << binary << " ";
		for(int j = 0; j < sample.size(); j++) {
			output_stream << sample.at(j) << " ";
		}
		output_stream << endl;
	}
}

void write_positional_encoding(vector<int> all_bp_positions, 
		vector<string> all_sample_IDs, 
		vector<vector<int>> smf, 
		string output_positional_encoding_file){
	// open output file to write encoding
        ofstream output_stream;
        output_stream.open(output_positional_encoding_file);

        // format sample ID float float float...
        // space delim
	int relative_position = -1;
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
		output_stream << all_sample_IDs[SID_i] << "_" << binary << " ";
                for(int j = 0; j < sample.size(); j++) {
			if (sample.at(j) == 1){
				relative_position = all_bp_positions.at(j) - all_bp_positions.at(0);
				output_stream << relative_position << " ";
			}

                }
                output_stream << endl;
        }

}
