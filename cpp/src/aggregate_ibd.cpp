#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <map>
#include <vector>

#include "aggregate_ibd.h"
#include "utils.h"

using namespace std;

int main(int argc, char* argv[]){
    string ibd_file = argv[1];
    char delim = ' ';
    main_ibd_aggregate(ibd_file, delim);

    
}

void main_ibd_aggregate(string ibd_file, char delim=' '){
/* open pairwise ibd file
 * 
 */
    ifstream ibd_file_stream;
    ibd_file_stream.open(ibd_file);
    cout << "Reading IBD file..." << ibd_file << endl;
    if (!ibd_file_stream.is_open()){
        cout << "FAILED TO OPEN: " << ibd_file << endl;
        exit(1);   
    }  

    string line;
    while (getline(ibd_file_stream, line)){
        vector<string> vec_line;
        split_line(line, delim, vec_line);
	cout << line[0] << endl;
    } 
    
}
