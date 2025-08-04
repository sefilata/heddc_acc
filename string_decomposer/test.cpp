#include <iostream>

#include "string_decomposer.h"
#include "../eddc_original/fasta.h"

using namespace std;


int main(int argc, char* argv[]){
	// TR配列
	// vector<int> s = {0,1,0,1,0,1};
	// vector<int> t = {0,1,0,1,0,1,0,1};
	vector<int> s = {'A','C','T','A','C','T','C','T','A','C','T'};
	// vector<int> s = {0,1,1,0,1,1,0,1,1,3,2,3,2,3,2,0,2,1};
	// vector<int> t = {0,1,1,0,1,1,3,2,3,2,3,2};
	// vector<int> s = {0,1,2,0,1,2,0,1,2};

	vector<vector<int>> units = {{'A','C','T'}};
	// vector<vector<int>> units = {{0,1,1}, {3,2}, {0,1}};
	// vector<vector<int>> units = {{0,1}};
	// vector<vector<int>> units = {{0,1,2}};

	vector<int> encoding;
	vector<vector<int>> decomposed_seq;
	SDParams params;

	string_decomposer(s, units, params, encoding, decomposed_seq);

	cout << "encoding: " << endl;
	for(int i : encoding){cout << i << " ";}
	cout << endl;
	cout << "decomposed_seq: " << endl;
	for(auto vec : decomposed_seq){
		cout << "(";
		for(int i = 0; i < vec.size(); i++){
			cout << vec[i];
			if(i < vec.size()-1){cout << " ";}
		}
		cout << ") ";
	}
	cout << endl;
}
