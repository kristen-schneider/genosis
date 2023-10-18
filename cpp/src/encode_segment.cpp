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

/**
 * encodes one vcf segment
 *
 * @param chrm_idx chromosome index
 * @param vcf_segment_file vcf segment file
 * @param sample_IDs_file file with all samples IDs to encode
 * @param encoding_def_map map whose key is a genotype and value is an encoding
 * @param chrm_bp_cm_map map whos key is a chromosome and whose value is a <bp, cm> map
 * @param output_gt_file output genotype encoding file
 * @param output_pos_file output positional encoding file
 *
 * @return void
 */
void encode_vectors(int chrm_idx,
		string vcf_segment_file,
		string sample_IDs_file,
                map<string, vector<int> > encoding_def_map, 
		map<int, map<int, float> > chrm_bp_cm_map,
		string output_gt_file,
		string output_pos_file){
	
    // converts vcf_segment_file from string to const char for htslib
    const char *vcf_segment_htslib = vcf_segment_file.c_str();
    
    // open VCF file with htslib
    htsFile *vcf_stream = bcf_open(vcf_segment_htslib, "r");
    if (!vcf_stream){
        cout << "FAILED TO OPEN: " << vcf_segment_htslib << endl;
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
    vector<vector<int> > all_gt_haplotype_encodings;
    // vector of positional haplotype encoding vectors
    vector<vector<int> > all_pos_haplotype_encodings;
    
    // vector of all bp in the VCF
    vector<int> all_bp_positions;
    // vector of all af in the VCF **not currently used**
    //vector<float> all_af;

    // READING GENOTYPES
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
        // **not currently used (below)
        //double_t qual = vcf_record->qual;				// quality
        //float *afs = 0;
        //int count = 0;
        //int ret = bcf_get_info_float(vcf_header, vcf_record, "AF", &afs, &count);
        //float af = afs[0];						// allele frequency
        
        all_bp_positions.push_back(pos+1);				// add pos to vector
        //all_af.push_back(af);						// add af to vector
        
        // reading genotypes
        string s_gt;        // string variable of genotype
        int *gt = NULL;     // int pointer to genotype
        int number_total_alleles = 0;
        int ngt_arr = 0;
        number_total_alleles = bcf_get_genotypes(vcf_header, vcf_record,  &gt, &ngt_arr);
        int alleles_per_gt = number_total_alleles/num_samples;
        
        // genotype encoding vector
        // < 0 0 0 1 0 1 1 1 1 > --> 
        // 5 samples with following gt encodings: 00 01 01 11 10
        vector<int> biallelic_gt_encoding_vector;
        
        // for each sample in the record, convert gt to encoding
        for (int i = 0; i < num_samples; i++){
        	
            // separate alleles for each gt
            int allele1 = bcf_gt_allele(gt[i*alleles_per_gt+0]);
            int allele2 = bcf_gt_allele(gt[i*alleles_per_gt+1]);
            // replace unknowns and concatinate genotypes to " | " format
            // if allele 1 is undefined by 0-3, unkown ('.')
            if (allele1 != 0 and allele1 != 1 and allele1 != 2 and allele1 != 3){
                // if allele 2 is undefined by 0-3, unkown ('.')
            	if (allele2 != 0 and allele2 != 1 and allele2 != 2 and allele2 != 3){
            	    s_gt = ".|.";    // both allels are unkown
            	}
            	else{
            	    if (allele2 < 0){
            		s_gt = ".|.";    // both alleles are unkown
            	    }
                    // if allele 2 is defined, convert string to int and store
            	    else{
            	        s_gt = ".|" + to_string(allele2);
            	    }
                }
            }
            // if allele 2 is undefined by 0-3, unkown ('.')
            if (allele2 != 0 and allele2 != 1 and allele2 != 2 and allele2 != 3){
                // if allele 1 is undefined by 0-3, unkown ('.')
            	if (allele1 != 0 and allele1 != 1 and allele1 != 2 and allele1 != 3){
            	    s_gt = ".|.";   // both alleles are unkown
            	}
            	else{
            	    if (allele1 < 0){
                        s_gt = ".|."; // both alleles are unkown
            	    }
            	    else{
            	    	s_gt = to_string(allele1) + "|.";
            	    }
                }
            }
                // both alleles are known and can be encoded 0-3
            else if (s_gt == ""){
            	s_gt = to_string(allele1)+"|"+to_string(allele2);
            }
            // mapping genotype to encoding vector and add to the final genotype vector
            // going from "0|1" --> 01
            vector<int> biallelic_encoding = encoding_def_map[s_gt];
            biallelic_gt_encoding_vector.push_back(biallelic_encoding[0]);
            biallelic_gt_encoding_vector.push_back(biallelic_encoding[1]);
            s_gt =  "";
        }
                
        // add current record to all genotype vectors	
        all_gt_haplotype_encodings.push_back(biallelic_gt_encoding_vector);
        biallelic_gt_encoding_vector.clear();

    } // end of reading records
    
    
    // transposing data frorm variant major to sample major format
    vector<vector<int> > sample_major_format_gt_vec = transpose_int(all_gt_haplotype_encodings); 
    
    // get list of sample IDs
    vector<string> all_sample_IDs = get_sample_IDs(sample_IDs_file); 
    
    // writing sample major format for genotype vectors
    write_SMF_gt(all_sample_IDs,
    		sample_major_format_gt_vec,
    		output_gt_file);
    
    // writing sample major format for positional vectors
    write_SMF_pos(chrm_idx,
    		all_sample_IDs,
    		sample_major_format_gt_vec,
    		all_bp_positions,
    		chrm_bp_cm_map,
    		output_pos_file);
}	

/**
 * write genotype sample major format out to file: sampleID 0 1 0 0 ...
 * 
 * @param all_sample_IDs all sample IDs for encoding vectors
 * @param smf sample major format for encoding vectors
 * @param output_gt_file output genotype encoding file
 *
 * @return void
 */
void write_SMF_gt(vector<string> all_sample_IDs,
		vector<vector<int> > smf,
		string output_gt_file){
	
    // open output file to write encoding
    ofstream output_stream;
    output_stream.open(output_gt_file);
    
    // format sample ID float float float...
    int SID_i = 0;
    int binary = -1;
    // for all sample haplotypes
    for (int i = 0; i < smf.size(); i++) {
        // at start, haploytype = 0
        if (binary == -1){
    	    binary = 0;
    	}
        // if haplotype 0 already covered, haplotype = 1
        else if (binary == 0){
    	    binary = 1;	// change flag to final 
    	}
        // if haplotyp 0 and 1 already covered, start for new sample
        else{
    	    SID_i ++;   // go to next sample ID
    	    binary = 0; // change to haplotype 0
    	}
        // current sample vector
    	vector<int> sample_vector = smf.at(i);
    	// write current sample ID with correct haplotype
        output_stream << all_sample_IDs[SID_i] << "_" << binary << " ";
        // for all records in this sample vector, write out
    	for(int j = 0; j < sample_vector.size(); j++) {
    	    output_stream << sample_vector.at(j) << " ";
    	}
    	output_stream << endl;
    }
}


/**
 * write positional sample major format out to file: sampleID 0 1 0 0 ...
 * 
 * @param chrm_idx chromosome index
 * @param all_sample_IDs all sample IDs for encoding vectors
 * @param smf_gt sample major format for genotype encoding vectors
 * @param all_bp_positions a lit of all basepair positions in a VCF segment
 * @param chrm_bp_cm_map map whos key is a chromosome and whose value is a <bp, cm> map
 * @param output_positional_encoding_file output positional encoding file
 *
 * @return void
 */
void write_SMF_pos(int chrm_idx,
		vector<string> all_sample_IDs, 
		vector<vector<int> > smf_gt, 
		vector<int> all_bp_positions,
		map<int, map<int, float> > chrm_bp_cm_map,
		string output_positional_encoding_file){
	
    // open output file to write encoding
    ofstream output_stream;
    output_stream.open(output_positional_encoding_file);
        
    // format sample ID float float float...
    // space delim
    int relative_position = -1;
    int SID_i = 0;
    int binary = -1;
    
    map<int, float> chromosome_chrm_bp_cm_map = chrm_bp_cm_map[chrm_idx];
    
    // for all sample haplotypes
    for (int i = 0; i < smf_gt.size(); i++) {
        // at start, haploytype = 0
        if (binary == -1){
    	    binary = 0;
    	}
        // if haplotype 0 already covered, haplotype = 1
        else if (binary == 0){
    	    binary = 1;	// change flag to final 
    	}
        // if haplotyp 0 and 1 already covered, start for new sample
        else{
    	    SID_i ++;   // go to next sample ID
    	    binary = 0; // change to haplotype 0
    	}
        // current sample vector
        vector<int> sample = smf_gt.at(i);
    	output_stream << all_sample_IDs[SID_i] << "_" << binary << " ";
    	
    	for(int j = 0; j < sample.size(); j++) {
            // only write positional encoding where there is a variant (1)
    	    if (sample.at(j) > 0){
    		int bp_pos = all_bp_positions.at(j);
    		float cm_pos = chromosome_chrm_bp_cm_map[bp_pos];
    		output_stream << cm_pos << " ";
    	    }
    
        }
        output_stream << endl;
    }
}

/* allele frequency not used for now
void write_SMF_af(vector<string> all_sample_IDs,
                vector<vector<int> > smf_gt,
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
				//float cm_pos = chrm_bp_cm_map[bp_pos];
				//output_stream << chrm_bp_cm_map[all_bp_positions.at(j)] << " ";
			}

                }
                output_stream << endl;
        }

}*/
