#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <math.h>
#include <map>
#include <string>
#include <vector>

#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <htslib/vcfutils.h>
#include <htslib/kstring.h>
#include <htslib/kseq.h>
#include <htslib/synced_bcf_reader.h>

#include "encode_positions.h"
#include "read_map.h"
#include "map_encodings.h"
#include "utils.h"


using namespace std;

void encode_positions(string map_file,
		string vcf_slice_file, 
		string sample_IDs_file,
		string output_encoding_file){
	
	// make map from map file
	map<int, float> bp_cm_map;
	bp_cm_map = make_bp_cm_map(map_file);
	/*
	for(auto it = bp_cm_map.cbegin(); it != bp_cm_map.cend(); ++it)
	{
    		std::cout << it->first << " " << it->second << "\n";
	}
	*/
	// converts vcfFile name to const char for htslib
	const char *vcf_slice = vcf_slice_file.c_str();
	// open VCF file with htslib
	cout << "VCF: " << vcf_slice << endl;
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

	vector<vector<float>> all_haplotype_pos_encodings;
	//vector<int> all_positions;
	
	cout << "...reading genotypes." << endl;
	while (bcf_read(vcf_stream, vcf_header, vcf_record) == 0){
		bcf_unpack(vcf_record, BCF_UN_ALL);
		bcf_unpack(vcf_record, BCF_UN_INFO);

		// making variable names for each column
		string chrm = bcf_hdr_id2name(vcf_header, vcf_record->rid);
		int pos = (unsigned long)vcf_record->pos + 1;
		string id = vcf_record->d.id;
		string ref = vcf_record->d.allele[0];
		string alt = vcf_record->d.allele[1];
		double_t qual = vcf_record->qual;
		//all_positions.push_back(pos);
		
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
		
		//map<string, vector<int>> encoding_map;
		//encoding_map = make_biallelic_encoding_map("/home/sdp/precision-medicine/example/example_encoding.txt");
		//vector<float> haplotype_0_pos_vector;
		//vector<float> haplotype_1_pos_vector;
		vector<float> haplotype_pos_encoding_vector;

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
			// record positions	
			//if (allele1 != 0 or allele2 != 0){
			//cout << allele1 << " " << allele2 << endl;
			float cm_pos = bp_cm_map[pos];
			//cout << pos << " " << cm_pos << endl;
			if (allele1 != 0){
				haplotype_pos_encoding_vector.push_back(cm_pos);
			}else if(allele1 == 0){
				haplotype_pos_encoding_vector.push_back(0);
			}
			if (allele2 != 0){
				haplotype_pos_encoding_vector.push_back(cm_pos);
			}else if(allele2 == 0){
                                haplotype_pos_encoding_vector.push_back(0);
                        }
			//}
			//vector<int> biallelic_encoding = encoding_map[s_gt];
			//cout << biallelic_encoding.size() << endl;
			//haplotype_pos_encoding_vector.push_back(biallelic_encoding[0]);
			//haplotype_pos_encoding_vector.push_back(biallelic_encoding[1]);
			// if non reference variant for first allele
			// if non reference variant for second allele
			//if (allele1 != 0 or allele2 != 0){

			//	float cm_pos = bp_cm_map[pos];		
				//cout << "allele2: " << allele2 << " " << pos << " " << cm_pos << endl;
			//	haplotype_pos_encoding_vector.push_back(cm_pos);
				//allele2 = bcf_gt_allele(gt[i*alleles_per_gt+1]);
				//s_gt = to_string(allele1)+"|"+to_string(allele2);
			// homozygous reference
			//else{
			//	continue;
			//}
			
			//vector<int> biallelic_encoding = encoding_map[s_gt];
			//haplotype_encoding_vector.push_back(biallelic_encoding[0]);
			//haplotype_encoding_vector.push_back(biallelic_encoding[1]);
				
		}
		
		//cout << haplotype_pos_encoding_vector.size() << endl;
		all_haplotype_pos_encodings.push_back(haplotype_pos_encoding_vector);
		/*cout << haplotype_pos_encoding_vector.size() << ": " << endl;
		for (int z = 0; z < haplotype_pos_encoding_vector.size(); z++){
			cout << haplotype_pos_encoding_vector[z] << " ";;
		}
		cout << endl;
		*/
		haplotype_pos_encoding_vector.clear();
		//cout << all_haplotype_pos_encodings.size() << endl;
	} // end of reading records
	cout << all_haplotype_pos_encodings.size() << endl;
	
	// transposing data
	cout << "...transposing data..." << endl;
	vector<vector<float>> sample_major_format_hap_pos_vec = transpose_float(all_haplotype_pos_encodings);

	vector<string> all_sample_IDs = get_sample_IDs(sample_IDs_file);
	
	// writing smf
	cout << "...writing sample major format encodings to file..." << endl;
	cout << sample_major_format_hap_pos_vec.size() << endl;
	write_SMF(all_sample_IDs, sample_major_format_hap_pos_vec, output_encoding_file);
}

/*
 * given a set of sample IDs,
 * and a set of smf encoded vectors, 
 * write: sampleID encoding
 */
void write_SMF(vector<string> all_sample_IDs,
	       	vector<vector<float>> smf,
	       	string output_encoding_file){
	
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
		
		vector<float> sample = smf.at(i);
		//cout << sample.size() << " ";
		//remove zeros
		//sample.erase(remove(sample.begin(), sample.end(), 0));
		//vector<float> sample_no_zeros;
		//sample_no_zeros = remove_zeros_float(sample);
		//cout << sample.size() << endl;
		// write to file
		output_stream << all_sample_IDs[SID_i] << "_" << binary << " ";
		for(int j = 0; j < sample.size(); j++) {
			float curr_encoding = sample.at(j);
			if(curr_encoding > 0){
				//cout << sample.at(j) << " ";
				output_stream << sample.at(j) << " ";
			}
		}
		//cout << endl;
		output_stream << endl;
	}
}

