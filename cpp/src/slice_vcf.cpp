#include <iostream>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>
#include <cstdlib>

#include "read_map.h"
#include "slice_vcf.h"
#include "utils.h"

using namespace std;

int slice_main(string map_file, int segment_size, string vcf_file, string out_base_name, string out_dir){
	// read map file and report number
	// of slices that should be generated
	vector<int> segment_SNP_counts;
	segment_SNP_counts = read_map_file(map_file, segment_size);
	int map_slice_count = segment_SNP_counts.size();
	cout << "...Counted " << map_slice_count << " " << segment_size << "cM slices.\n" << endl;

	// read header of full chromosome VCF file
	// store as list to write to header of 
	// smaller VCF files
	cout << "...Reading VCF file..." << endl;
	vector<string> vcf_header = read_vcf_header(vcf_file);
	int header_num_lines = vcf_header.size();
	cout << "...Read " << header_num_lines << " lines from VCF header." << endl;


	// slice full chromosome VCF file into smaller slices
	cout << "...Slice size: " << segment_size << "cM." << endl;
	int num_segments = slice(vcf_file, vcf_header, segment_SNP_counts,
        	out_base_name, out_dir);
	
	cout << "...Done reading VCF file." << endl;
	cout << "Done slicing." << endl;
	return num_segments;
}

/*
 * Open full chromosome VCF file and read the header.
 * Store the header data in a vector of strings
 * Return the header data. (To be written at the
 * top of every smaller VCF file)
 */
vector<string> read_vcf_header(string vcf_file){
	// stores VCF header
	vector<string> vcf_header;
	
	// open file and check success
    	ifstream vcf_file_stream;
    	vcf_file_stream.open(vcf_file);
    	if (!vcf_file_stream.is_open()){
        	cout << "FAILED TO OPEN: " << vcf_file << endl;
		exit(1);
    	}
        
	// read vcf file header and stop
	// when header is done
	string line;
        while (getline (vcf_file_stream, line)){
            char char1 = line.at(0);
            if (char1 == '#'){
                vcf_header.push_back(line);
            }
	    // stop when header is done
            else{ break; }
        }	
	return vcf_header;
}

/*
 * Open full VCF file, ignore header, 
 * write a one slice to slice file
 * return number of slices
 */
int slice(string vcf_file, 
		vector<string> vcf_header, 
		vector<int> segment_SNP_counts, 
		string base_name, string out_dir){

	cout << out_dir << endl;
	// to return
	int slice_index = 0;
	int slice_snp_count = segment_SNP_counts[slice_index];
	//int slice_end = segment_SNP_counts[slice_index];
    	// open vcf file and check success
    	ifstream vcf_file_stream;
    	vcf_file_stream.open(vcf_file);
    	if (!vcf_file_stream.is_open()){
    		cout << "FAILED TO OPEN: " << vcf_file << endl;
        	exit(1);
	}
        int total_line_count = 0;
        int SNPS_in_slice = 0;
	
	// at start of a new slice,
	// open new slice file
	// write header
	ofstream slice_file_stream;
	// name slice out file
        string out_vcf_slice_file = out_dir + base_name + \
        	".seg." + to_string(slice_index) + \
                ".vcf";
	//cout << "... ... slice " << slice_count << endl;
       	slice_file_stream.open(out_vcf_slice_file);
        for (int i = 0; i < vcf_header.size(); i ++){
                slice_file_stream << vcf_header[i] << endl;
        }

	// read vcf file until end
	string line;
        while (getline (vcf_file_stream, line)){
            	/*
		// FOR TESTING
		if (slice_count >= 10){
			break;
		}
		*/
		// ignore header
		char char1 = line.at(0);
            	if (char1 == '#'){continue;}
		else{
			// still building a slice
			if (SNPS_in_slice < slice_snp_count){
				slice_file_stream << line << endl;
				SNPS_in_slice ++;
				total_line_count ++;
			}
			// reached end of slice-->increment slice count
			// close file-->increment slice count
			// open new file-->write header-->write line
			else if (SNPS_in_slice == slice_snp_count){
				//slice_file_stream << line;
				slice_file_stream.close();
				slice_index += 1;

				// open next slice file and write header 
                                string out_vcf_slice_file = out_dir + base_name + \
                                        ".seg." + to_string(slice_index) + \
                                        ".vcf";
                                slice_file_stream.open(out_vcf_slice_file);
                                for (int i = 0; i < vcf_header.size(); i ++){
                                        slice_file_stream << vcf_header[i] << endl;
                                }
				slice_file_stream << line << endl;
				SNPS_in_slice = 1;
				slice_snp_count = segment_SNP_counts[slice_index];
			}
		}
        }
	// write last line
	if (SNPS_in_slice < slice_snp_count && line.size() > 0){
		slice_file_stream << line << endl;
	}
	slice_file_stream.close();	
        
	return slice_index;
}
/*
 * Open a map file and create a vector of 
 * start basepairs for each segment
 
vector<int> segment_end_bp(string map_file, float slice_size){
	int max_cm = slice_size; // start max at slice sice
	int cm_index = 2;	 // column index where cm data is
	int bp_index = 3;	 // column index where bp data is 
	vector<int> slice_end_points;	// vector of all ending points

	ifstream map_file_stream;
	map_file_stream.open(map_file);
	if (!map_file_stream.is_open()){
		cout << "FAILED TO OPEN: " << map_file << endl;
		exit(1);
	}
	cout << "...Reading map file..." << endl;
	string line;
	while (getline (map_file_stream, line)){
		split_line(line, ' ', single_SNP);
		float record_cm = stof(single_SNP[cm_index]);
		int record_bp = stoi(single_SNP[bp_index]);
		if (record_cm >= max_cm){
			slice_end_points.push_back(record_bp);
			max_cm += slice_size;
		}
	}
	return slice_end_points;
}
*/

