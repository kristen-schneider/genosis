#include <cstdlib>
#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <vector>

#include "encode_positions.h"
#include "read_map.h"


using namespace std;

void encode_positions(string map_file,
		string vcf_slice_file){
	// make map from map file
	map<int, float> bp_cm_map;
	bp_cm_map = make_bp_cm_map(map_file);

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

	vector<vector<int>> all_haplotype_encodings;
	vector<int> all_positions;
	
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
		all_positions.push_back(pos);
		
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
		vector<int> haplotype_encoding_vector;

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
			vector<int> biallelic_encoding = encoding_map[s_gt];
			haplotype_encoding_vector.push_back(biallelic_encoding[0]);
			haplotype_encoding_vector.push_back(biallelic_encoding[1]);
				
		}
		all_haplotype_encodings.push_back(haplotype_encoding_vector);

		haplotype_encoding_vector.clear();

}
