#include <iostream>
#include <map>

using namespace std;

map<string,int> readEncodingFile(string encodingFile){
	map<string, int> encodingMap;
	return encodingMap;
}

int validateEncodingCount(map<string,int> m){
	return m.size();
}

int main(void){
	string encodingFile = "/home/sdp/precision-medicine/encodings/trivial.txt";
	map<string, int> m = readEncodingFile(encodingFile);
	validateEncodingCount(m);
}
