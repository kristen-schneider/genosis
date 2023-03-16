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

#include "encode_segment.h"
#include "map_encodings.h"
#include "read_config.h"
#include "utils.h"

using namespace std;

/*
 * takes one sliced vcf file and encoding instructions
 * and writes output encoiding file
 */
void encode_vectors(int chrm_idx,
		string vcf_slice_file,
		string sample_IDs_file, 
		map<string, vector<int>> encoding_map, 
		map<int, map<int, float>> bp_cm_map,
		string output_gt_file,
		string output_pos_file, 
		string output_af_file){
	

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

	// vector of nomical haplotype encoding vectors
	vector<vector<int>> all_gt_haplotype_encodings;
	// vector of positional haplotype encoding vectors
	vector<vector<float>> all_pos_haplotype_encodings;
	
	// vector of all bp in the VCF
	vector<int> all_bp_positions;
	// vector of all af in the VCF
	vector<float> all_af;

	cout << "...Reading genotypes." << endl;
	while (bcf_read(vcf_stream, vcf_header, vcf_record) == 0){
		bcf_unpack(vcf_record, BCF_UN_ALL);
		bcf_unpack(vcf_record, BCF_UN_INFO);

		// storing data for each column
		string chrm = bcf_hdr_id2name(vcf_header, vcf_record->rid);	// chrom (string)
		int chrm_int = stoi(chrm.substr(3, chrm.length()));		// chrom (int)
		int pos = vcf_record->pos;
		//int pos = (unsigned long)vcf_record->pos;			// pos
		string id = vcf_record->d.id;					// sample id
		string ref = vcf_record->d.allele[0];				// ref allele
		string alt = vcf_record->d.allele[1];				// alt allele
		double_t qual = vcf_record->qual;				// quality
		float *afs = 0;
		int count = 0;
		int ret = bcf_get_info_float(vcf_header, vcf_record, "AF", &afs, &count);
		float af = afs[0];						// allele frequency
	
		all_bp_positions.push_back(pos+1);				// add pos to vector
		all_af.push_back(af);						// add af to vector

		// reading genotypes
		string s_gt;
		int *gt = NULL;
		int number_total_alleles = 0;
		int ngt_arr = 0;
		number_total_alleles = bcf_get_genotypes(vcf_header, vcf_record,  &gt, &ngt_arr);
		int alleles_per_gt = number_total_alleles/num_samples;

		// single variant vector
		vector<int> haplotype_gt_encoding_vector;
		vector<int> haplotype_pos_encoding_vector;
		
		// for each sample in the record, convert gt to encoding
		for (int i = 0; i < num_samples; i++){
			
			// separate two alleles for each gt
			int allele1 = bcf_gt_allele(gt[i*alleles_per_gt+0]);
			int allele2 = bcf_gt_allele(gt[i*alleles_per_gt+1]);
			// replace unknowns and concatinate genotypes to " | " format
			if (allele1 != 0 and allele1 != 1 and allele1 != 2 and allele1 != 3){
				if (allele2 != 0 and allele2 != 1 and allele2 != 2 and allele2 != 3){
					s_gt = ".|.";
				}
				else{
					s_gt = ".|" + to_string(allele2);
				}
			}
			if (allele2 != 0 and allele2 != 1 and allele2 != 2 and allele2 != 3){
				if (allele1 != 0 and allele1 != 1 and allele1 != 2 and allele1 != 3){
					s_gt = ".|.";
				}
				else{
					s_gt = to_string(allele1) + "|.";
				}
                        }
			else{
				s_gt = to_string(allele1)+"|"+to_string(allele2);
			}
			// mapping genotype to encoding vector and add to the hapltype vector
			vector<int> biallelic_encoding = encoding_map[s_gt];
			haplotype_gt_encoding_vector.push_back(biallelic_encoding[0]);
			haplotype_gt_encoding_vector.push_back(biallelic_encoding[1]);
			// record the cm position for the current record (variant)
			// add cm value if variant allele, add 0 if ref allele
		}
		
		all_gt_haplotype_encodings.push_back(haplotype_gt_encoding_vector);
		all_pos_haplotype_encodings.push_back(haplotype_pos_encoding_vector);
		haplotype_gt_encoding_vector.clear();

	} // end of reading records
	
	cout << "...Done reading genotypes." << endl;

	// transposing data
	cout << "...Transposing data." << endl;
	//cout << "......gt data......" << endl;
	vector<vector<int>> sample_major_format_gt_vec = transpose_int(all_gt_haplotype_encodings); 
	//cout << "......pos data......" << endl;
	//vector<vector<float>> sample_major_format_pos_vec = transpose_float(all_pos_haplotype_encodings);
	cout << "...Done transposing data." << endl;

	vector<string> all_sample_IDs = get_sample_IDs(sample_IDs_file); 
	// writing smf
	cout << "...Writing gt encodings to file: " << output_gt_file << "." << endl;
	write_SMF_gt(all_sample_IDs,
			sample_major_format_gt_vec,
			output_gt_file);

	cout << "...Writing pos encodings to file: " << output_pos_file << "." << endl;
	write_SMF_pos(chrm_idx,
			all_sample_IDs,
			sample_major_format_gt_vec,
			all_bp_positions,
			bp_cm_map,
			output_pos_file);
	/*
	cout << "...Writing af encodings to file: " << output_af_file << "." << endl;
	write_SMF_af(all_sample_IDs,
			sample_major_format_gt_vec,
			all_af,
			output_af_file);
	*/
	cout << "...Done writing sample major format." << endl;

	// writing positinal encoding
	//cout << "...writing positional encodings to file..." << endl;
	//write_positional_encoding(all_bp_positions, all_sample_IDs, sample_major_format_gt_vec, output_position_file);	
}	

/*
 * given a set of sample IDs,
 * and a set of smf encoded vectors, 
 * write: sampleID encoding
 */
void write_SMF_gt(vector<string> all_sample_IDs,
		vector<vector<int>> smf,
		string output_gt_file){
	
	// open output file to write encoding
	ofstream output_stream;
	output_stream.open(output_gt_file);

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

void write_SMF_pos(int chrm_idx,
		vector<string> all_sample_IDs, 
		vector<vector<int>> smf_gt, 
		vector<int> all_bp_positions,
		map<int, map<float, float>> bp_cm_map,
		string output_positional_encoding_file){
	
	// open output file to write encoding
        ofstream output_stream;
        output_stream.open(output_positional_encoding_file);

	
        // format sample ID float float float...
        // space delim
	int relative_position = -1;
        int SID_i = 0;
        int binary = -1;

	map<float, float> chromosome_bp_cm_map[chrm_idx];

        //cout << "bp_num: " << all_bp_positions.size() << endl;
	for (int i = 0; i < smf_gt.size(); i++) {
                if (binary == -1){
                        binary = 0;
                }else if (binary == 0){
                        binary = 1;
                }else{
                        SID_i ++;
                        binary = 0;
                }
                vector<int> sample = smf_gt.at(i);
		output_stream << all_sample_IDs[SID_i] << "_" << binary << " ";
		
		for(int j = 0; j < sample.size(); j++) {
			if (sample.at(j) > 0){
				float bp_pos = all_bp_positions.at(j);
				float cm_pos = chromosome_bp_cm_map[bp_pos];
				//float cm_pos = bp_cm_map[bp_pos];
				output_stream << bp_cm_map[all_bp_positions.at(j)] << " ";
			}

                }
                output_stream << endl;
        }
	

}

void write_SMF_af(vector<string> all_sample_IDs,
                vector<vector<int>> smf_gt,
                vector<float> all_af,
                string output_af_file){
	// open output file to write encoding
        ofstream output_stream;
        output_stream.open(output_af_file);

        // format sample ID float float float...
        // space delim
        int SID_i = 0;
        int binary = -1;
	for (int i = 0; i < smf_gt.size(); i++) {
                if (binary == -1){
                        binary = 0;
                }else if (binary == 0){
                        binary = 1;
                }else{
                        SID_i ++;
                        binary = 0;
                }
                vector<int> sample = smf_gt.at(i);
		output_stream << all_sample_IDs[SID_i] << "_" << binary << " ";
		
		for(int j = 0; j < sample.size(); j++) {
			if (sample.at(j) > 0){
				float af = all_af.at(j);
				output_stream << af << " ";
				//float bp_pos = all_bp_positions.at(j);
				//float cm_pos = bp_cm_map[bp_pos];
				//output_stream << bp_cm_map[all_bp_positions.at(j)] << " ";
			}

                }
                output_stream << endl;
        }

}
