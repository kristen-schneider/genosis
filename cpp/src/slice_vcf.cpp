#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "read_config.h"
#include "slice_vcf.h"

using namespace std;

int main(int argc, char* argv[]){

	// read config file for chromosome VCF file
	string configFile = argv[1];   // configuration file will all options
	map<string, string> config_options;
	config_options = get_config_options(configFile);

	string vcf_file = config_options["vcf_file"];
	int slice_size = stoi(config_options["slice_size"]);
	string out_dir = config_options["out_dir"];
    	string out_base_name = config_options["out_base_name"];

	// read header of full chromosome VCF file
	// store as list to write to header of 
	// smaller VCF files
	vector<string> vcf_header = read_vcf_header(vcf_file);
	int header_num_lines = vcf_header.size();
	cout << "Read " << header_num_lines << " lines from VCF header." << endl;


	// slice full chromosome VCF file into smaller slices
	int num_slices = slice(vcf_file, vcf_header, slice_size,
        	out_base_name, out_dir);
	return 0;
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
int slice(string vcf_file, vector<string> vcf_header, 
		int slice_size, 
		string base_name, string out_dir){

	// to return
	int slice_count = 0;

    	// open vcf file and check success
    	ifstream vcf_file_stream;
    	vcf_file_stream.open(vcf_file);
    	if (!vcf_file_stream.is_open()){
    		cout << "FAILED TO OPEN: " << vcf_file << endl;
        	exit(1);
	}
        int total_line_count = 0;
        int lines_in_slice = 0;
	
	// at start of a new slice,
	// open new slice file
	// write header
	ofstream slice_file_stream;
	// name slice out file
        string out_vcf_slice_file = out_dir + base_name + \
        	".seg." + to_string(slice_count) + \
                ".vcf";
       	slice_file_stream.open(out_vcf_slice_file);
        for (int i = 0; i < vcf_header.size(); i ++){
                slice_file_stream << vcf_header[i] << endl;
        }

	// read vcf file until end
	string line;
        while (getline (vcf_file_stream, line)){
            	// ignore header
		char char1 = line.at(0);
            	if (char1 == '#'){continue;}
		else{
			// still building a slice
			if (lines_in_slice < slice_size){
				slice_file_stream << line << endl;
				lines_in_slice ++;
				total_line_count ++;
			}
			// reached end of slice-->increment slice count
			// close file-->increment slice count
			// open new file-->write header-->write line
			else if (lines_in_slice == slice_size){
				//slice_file_stream << line;
				slice_file_stream.close();
				slice_count += 1;
				
				// open next slice file and write header 
                                string out_vcf_slice_file = out_dir + base_name + \
                                        ".seg." + to_string(slice_count) + \
                                        ".vcf";
                                slice_file_stream.open(out_vcf_slice_file);
                                for (int i = 0; i < vcf_header.size(); i ++){
                                        slice_file_stream << vcf_header[i] << endl;
                                }
				slice_file_stream << line << endl;
				lines_in_slice = 1;	
			}
		}
        }
	// write last line
	if (lines_in_slice < slice_size){
		slice_file_stream << line << endl;
	}
	slice_file_stream.close();
		
        slice_file_stream.close();
        cout << "SEGMENT LENGTH: " << slice_size << endl;
        cout << "LINE COUNT: " << total_line_count << endl;
        cout << "NUM SEGMENTS: " << slice_count << endl;
        int num_variants = total_line_count;
	return slice_count;
}
