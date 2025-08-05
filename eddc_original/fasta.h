#pragma once
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <map>
#include <set>

using namespace std;


// [A, T, G, C] -> [0, 1, 2, 3]
map<char, int> base_map = {
	{'A', 0}, {'T', 1}, {'G', 2}, {'C', 3},
	{'a', 0}, {'t', 1}, {'g', 2}, {'c', 3},
	{'0', 0}, {'1', 1}, {'2', 2}, {'3', 3}	// これはデバッグ用
};
map<int, char> base_map_reverse = {
	{0, 'A'}, {1, 'T'}, {2, 'G'}, {3, 'C'}
};
void string_to_digits(const string &line, vector<int> &read){
	for(int i = 0; i < (int)line.length(); i++){
		auto it = base_map.find(line[i]);
		if(it == base_map.end()){
			cerr << "[Error] (string_to_digits): no matching base character: " << line[i] << endl;
		}
		read[i] = it->second;
	}
}
// 関数名おかしくない？
void encoding_to_digits(const string &line, vector<int> &encoding){
	encoding.clear();
	istringstream iss(line);
	int v;
	while (iss >> v) {
		encoding.push_back(v);
	}
}
void digits_to_string(const vector<int> &read, string &line){
	line.clear();
	for(int i = 0; i < (int)read.size(); i++){
		auto it = base_map_reverse.find(read[i]);
		if(it == base_map_reverse.end()){
			cerr << "[Error] (digits_to_string): no matching base value: " << read[i] << endl;
		}
		line += it->second;
	}
}

// 2025.08.14 同一ユニットの除去を実装
void read_fasta(
	const string &filepath,
	vector<string> &read_names,
	vector<vector<int>> &reads,
	vector<vector<int>> &units
){
	ifstream infile(filepath);
	if(!infile.is_open()){
		cerr << "Error opening file: " << filepath << endl;
		return;
	}

	// 次が read 行か unit 行か（read -> 'r', unit -> 'u'）
	// read 行: "> #", unit 行: "> f"
	// 各リード（ユニット）は配列内で改行されていない
	char next = 'r';

	string line;
	set<string> units_set;
	while(getline(infile, line)){
		if(line.empty()){
			continue;
		}

		if(line[0] == '>'){
			if(line[1] == ' ' && line[2] == 'f'){
				next = 'u';
			}else if(line[1] == ' ' && line[2] == '#'){
				next = 'r';
				read_names.push_back(line);
			}else{
				cerr << "[Error] : unexpected header in input fasta" << endl;
			}
		}else{
			vector<int> read(line.length());
			string_to_digits(line, read);
			if(next == 'r'){reads.push_back(read);}
			else{
				if(units_set.count(line) == 0){
					units_set.insert(line);
					units.push_back(read);
				}
			}
		}
	}

	infile.close();
}

// 2025.01.02	リードファイルとユニットファイルを分離
// 2025.05.22	ユニットの重複を除去する操作を追加
void read_fasta2(
	const string &readfile,
	const string &unitfile,
	vector<string> &read_names,
	vector<vector<int>> &reads,
	vector<vector<int>> &units
){
	// read file
	ifstream infile(readfile);
	if(!infile.is_open()){
		cerr << "Error opening read file: " << readfile << endl;
		return;
	}
	// 各リード（ユニット）は配列内で改行されていない
	string line;
	while(getline(infile, line)){
		if(line.empty()){
			continue;
		}

		if(line[0] == '>'){
			read_names.push_back(line);
		}else{
			vector<int> read(line.length());
			string_to_digits(line, read);
			reads.push_back(read);
		}
	}
	infile.close();

	// unit file
	set<string> units_set;
	ifstream infile2(unitfile);
	if(!infile2.is_open()){
		cerr << "Error opening unit file: " << unitfile << endl;
		return;
	}
	// 各リード（ユニット）は配列内で改行されていない
	while(getline(infile2, line)){
		if(line.empty()){
			continue;
		}

		if(line[0] != '>'){
			if (units_set.count(line) == 0) {
				units_set.insert(line);
				vector<int> unit(line.length());
				string_to_digits(line, unit);
				units.push_back(unit);
			}
		}
	}
	infile2.close();
}

// 2025.05.02	encodingfile(半角スペース区切りユニット番号列)を追加で受け取る
void read_fasta3(
	const string &readfile,
	const string &unitfile,
	const string &encodingfile,
	vector<string> &read_names,
	vector<vector<int>> &reads,
	vector<vector<int>> &units,
	vector<vector<int>> &encodings
){
	// read file
	ifstream infile(readfile);
	if(!infile.is_open()){
		cerr << "Error opening read file: " << readfile << endl;
		return;
	}
	// 各リード（ユニット）は配列内で改行されていない
	string line;
	while(getline(infile, line)){
		if(line.empty()){
			continue;
		}

		if(line[0] == '>'){
			read_names.push_back(line);
		}else{
			vector<int> read(line.length());
			string_to_digits(line, read);
			reads.push_back(read);
		}
	}
	infile.close();

	// encoding file
	ifstream infile2(encodingfile);
	if(!infile2.is_open()){
		cerr << "Error opening encoding file: " << encodingfile << endl;
		return;
	}
	while(getline(infile2, line)){
		if(line.empty()){
			continue;
		}

		if(line[0] != '>'){
			vector<int> encoding;
			encoding_to_digits(line, encoding);
			encodings.push_back(encoding);
		}
	}
	infile2.close();

	// unit file
	ifstream infile3(unitfile);
	if(!infile3.is_open()){
		cerr << "Error opening unit file: " << unitfile << endl;
		return;
	}
	// 各リード（ユニット）は配列内で改行されていない
	while(getline(infile3, line)){
		if(line.empty()){
			continue;
		}

		if(line[0] != '>'){
			vector<int> unit(line.length());
			string_to_digits(line, unit);
			units.push_back(unit);
		}
	}
	infile3.close();
}
