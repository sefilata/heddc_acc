#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <math.h>
#include <chrono>
#include <cstdlib>
#include <climits>

#include "../eddc_original/fasta.h"
#include "hEDDC_par1.h"
#include "../string_decomposer/string_decomposer.h"

using namespace std;

// Pattern of parameters of dup/cont; 0: 0.5*unit_length, 1: all 0.5
const int PARAMETER_PATTERN = 0;


void string_decompose(
	const vector<vector<int>> &reads,
	const vector<string> &read_names,
	const vector<vector<int>> &units,
	vector<vector<int>> &encodings,
	vector<vector<vector<int>>> &decomposed_seqs,
	const string &encodings_file
){
	int tr_num = reads.size();
	encodings.resize(tr_num);
	decomposed_seqs.resize(tr_num);

	SDParams sdparams;
	for(int i = 0; i < tr_num; i++){
		string_decomposer(reads[i], units, sdparams, encodings[i], decomposed_seqs[i]);
	}

	// encodingsの出力 (既にファイルが存在する場合はスキップ)
	if(!encodings_file.empty()){
		ofstream ofs(encodings_file);
		if(!ofs.is_open()){
			cerr << "Error opening encodings file: " << encodings_file << endl;
			return;
		}
		ofs << "# units: ";
		for(int i = 0; i < units.size(); i++){
			ofs << "(" << i << ", ";
			string unit_str;
			digits_to_string(units[i], unit_str);
			ofs << unit_str;
			ofs << ") ";
		}
		ofs << endl;
		for(int i = 0; i < tr_num; i++){
			ofs << read_names[i] << endl;
			for(int j = 0; j < encodings[i].size(); j++){
				ofs << encodings[i][j] << " ";
			}
			ofs << endl;
		}
		ofs.close();
	}
}

// スコアをtsv形式で標準出力
void out_scores(
	const vector<vector<Score>> &scores,
	const vector<vector<int>> &reads,
	string out_file
){
	int read_num = reads.size();

	ofstream ofs(out_file);
	if(!ofs.is_open()){
		cerr << "Error opening output file: " << out_file << endl;
		return;
	}
	for(int i = 0; i < read_num; i++){
		for(int j = 0; j < read_num; j++){
			ofs << setprecision(8) << scores[i][j].get_score() / sqrt(reads[i].size() * reads[j].size());;
			if(j != read_num-1) ofs << "\t";
		}
		ofs << endl;
	}
}

// 変異情報の出力
void out_variants(
	const vector<vector<Score>> &scores,
	const string &out_file,
	const vector<vector<int>> &units
){
	ofstream ofs(out_file);
	if(!ofs.is_open()){
		cerr << "Error opening variant output file: " << out_file << endl;
		return;
	}

	for(int i = 0; i < scores.size(); i++){
		for(int j = 0; j < scores[i].size(); j++){
			ofs << "{mut:" << scores[i][j].get_mut();
			ofs << ", indel:" << scores[i][j].get_indel();
			ofs << ", dup:";
			scores[i][j].print_dup(ofs, units);
			ofs << "}";
			if(j != scores[i].size()-1) ofs << "\t";
		}
		ofs << endl;
	}
	ofs.close();
}

// 実行時間の出力
void out_time(long long dur_msec, long long dur_usec, string out_file, const vector<long long> &measure_time){
	ofstream ofs(out_file);
	if(!ofs.is_open()){
		cerr << "Error opening time output file: " << out_file << endl;
		return;
	}
	ofs << dur_msec << " msec" << endl;
	ofs << dur_usec << " usec" << endl;
	ofs << measure_time[0] << " msec (c1,c2,valid_rules)" << endl;
	ofs << measure_time[1] << " msec (f_scores)" << endl;
	ofs << measure_time[2] << " msec (main dp)" << endl;
	ofs.close();
}

int main(int argc, char* argv[]){
	if(argc < 6){
		cerr << "Usage: " << argv[0] << " <read_unit_fasta> <score_tsv> <variant_tsv> <time_tsv> <encodings_fasta>" << endl;
		return 1;
	}
	string read_unit_file = argv[1]; 	// uTRの出力fasta
	string score_file = argv[2]; 		// EDDCのスコア出力ファイル
	string variant_file = argv[3]; 		// 変異情報の出力ファイル
	string time_file = argv[4]; 		// 実行時間の出力ファイル
	string encodings_file = argv[5]; 	// string decomposerの出力ファイル (既にある場合は出力しない)

	// fastaのパース, 同一ユニットの除去あり
	vector<string> read_names;
	vector<vector<int>> reads;
	vector<vector<int>> units;
	read_fasta(read_unit_file, read_names, reads, units);

	// mut, indel, dupのスコア（eddc_units は共通のものを使えるようにあとで修正）
	vector<vector<int>> eddc_units = {{0}, {1}, {2}, {3}};
	for(auto vec : units){eddc_units.push_back(vec);}
	Params params(1.0, 1.0, PARAMETER_PATTERN, eddc_units);

	// string decomposer
	vector<vector<int>> encodings;
	vector<vector<vector<int>>> decomposed_seqs;
	string_decompose(reads, read_names, units, encodings, decomposed_seqs, encodings_file);

	// EDDCの計算
	vector<vector<Score>> scores;
	vector<long long> measure_time;
	auto start = chrono::high_resolution_clock::now();
	heddc_all(encodings, units, params, scores, measure_time);
	auto end = chrono::high_resolution_clock::now();
	auto dur_msec = chrono::duration_cast<chrono::milliseconds>(end - start).count();
	auto dur_usec = chrono::duration_cast<chrono::microseconds>(end - start).count();

	// 出力
	out_scores(scores, reads, score_file);
	out_variants(scores, variant_file, units);
	out_time(dur_msec, dur_usec, time_file, measure_time);
}