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
extern map<char, int> base_map;
extern map<int, char> base_map_reverse;
void string_to_digits(const string &line, vector<int> &read);
// 関数名おかしくない？
void encoding_to_digits(const string &line, vector<int> &encoding);
void digits_to_string(const vector<int> &read, string &line);

// 2025.08.14 同一ユニットの除去を実装
void read_fasta(
	const string &filepath,
	vector<string> &read_names,
	vector<vector<int>> &reads,
	vector<vector<int>> &units
);

// 2025.01.02	リードファイルとユニットファイルを分離
// 2025.05.22	ユニットの重複を除去する操作を追加
void read_fasta2(
	const string &readfile,
	const string &unitfile,
	vector<string> &read_names,
	vector<vector<int>> &reads,
	vector<vector<int>> &units
);

// 2025.05.02	encodingfile(半角スペース区切りユニット番号列)を追加で受け取る
void read_fasta3(
	const string &readfile,
	const string &unitfile,
	const string &encodingfile,
	vector<string> &read_names,
	vector<vector<int>> &reads,
	vector<vector<int>> &units,
	vector<vector<int>> &encodings
);
