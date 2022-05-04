#include <iostream>
#include <fstream>
#include <vector>
#include "utils.h"
#include "slice.h"


using namespace std;

int main(void){

	int x = 0;
	cout << x << endl;
	cout << "Start of program." << endl;
	cout << "...reading vcf..." << endl;
	sliceVCF();
	cout << "End of program." << endl;
	return 0;
}
