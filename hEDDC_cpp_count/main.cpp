#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <math.h>
#include <chrono>
#include <cstdlib>
#include <climits>

#include "../eddc_original/fasta.h"
#include "hEDDC_acc_par1.h"
#include "../string_decomposer/string_decomposer.h"

using namespace std;

int SCORE_PRECISION = 8; // Number of decimal places for score output


// String decomposer and encodings output
void string_decompose(
	const vector<vector<int>> &reads,
	const vector<string> &read_names,
	const vector<vector<int>> &units,
	vector<vector<int>> &encodings,
	vector<vector<vector<int>>> &decomposed_seqs,
	const string &encodings_file,
	int e_flag 	// Whether to output encodings file or not
){
	int tr_num = reads.size();
	encodings.resize(tr_num);
	decomposed_seqs.resize(tr_num);

	SDParams sdparams;
	for(int i = 0; i < tr_num; i++){
		string_decomposer(reads[i], units, sdparams, encodings[i], decomposed_seqs[i]);
	}

	// encodingsの出力
	if(e_flag){
		ofstream ofs(encodings_file);
		if(!ofs.is_open()){
			cerr << "Error opening encodings file: " << encodings_file << endl;
			return;
		}
		ofs << "# units: ";
		for(int i = 0; i < (int)units.size(); i++){
			ofs << "(" << i << ", ";
			string unit_str;
			digits_to_string(units[i], unit_str);
			ofs << unit_str;
			ofs << ") ";
		}
		ofs << endl;
		for(int i = 0; i < tr_num; i++){
			ofs << read_names[i] << endl;
			for(int j = 0; j < (int)encodings[i].size(); j++){
				ofs << encodings[i][j] << " ";
			}
			ofs << endl;
		}
		ofs.close();
	}
}

// Score output in tsv file
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
			ofs << setprecision(SCORE_PRECISION) << scores[i][j].get_score() / sqrt(reads[i].size() * reads[j].size());;
			if(j != read_num-1) ofs << "\t";
		}
		ofs << endl;
	}
}

// Score output in stdout
void print_scores(
	const vector<vector<Score>> &scores,
	const vector<vector<int>> &reads
){
	int read_num = reads.size();

	for(int i = 0; i < read_num; i++){
		for(int j = 0; j < read_num; j++){
			cout << setprecision(SCORE_PRECISION) << scores[i][j].get_score() / sqrt(reads[i].size() * reads[j].size());;
			if(j != read_num-1) cout << "\t";
		}
		cout << endl;
	}
}

// Output variations (mutations, indels, duplications/contractions)
void out_variants(
	const vector<vector<Score>> &scores,
	const string &out_file
){
	ofstream ofs(out_file);
	if(!ofs.is_open()){
		cerr << "Error opening variant output file: " << out_file << endl;
		return;
	}

	for(int i = 0; i < (int)scores.size(); i++){
		for(int j = 0; j < (int)scores[i].size(); j++){
			ofs << "{mut:" << scores[i][j].get_mut();
			ofs << ", indel:" << scores[i][j].get_indel();
			ofs << ", dup:";
			scores[i][j].print_dup(ofs);
			ofs << "}";
			if(j != (int)scores[i].size()-1) ofs << "\t";
		}
		ofs << endl;
	}
	ofs.close();
}

// Output execution time
void out_time(long long dur_msec, long long dur_usec, string out_file, const vector<long long> &measure_time){
	ofstream ofs(out_file);
	if(!ofs.is_open()){
		cerr << "Error opening time output file: " << out_file << endl;
		return;
	}
	ofs << dur_msec << " msec" << endl;
	// ofs << dur_usec << " usec" << endl;
	ofs << measure_time[0] << " msec (c1,c2,valid_rules)" << endl;
	ofs << measure_time[1] << " msec (f_scores)" << endl;
	ofs << measure_time[2] << " msec (main dp)" << endl;
	ofs.close();
}

int main(int argc, char* argv[]){
	string read_file; 		// TR fasta (in)
	string unit_file; 		// Unit fasta (in)
	string score_file; 		// Score file (out)
	string variant_file; 	// Variations file (out)
	string time_file; 		// Execution time file (out)
	string encodings_file; 	// Etring decomposition results file (out)

	// Parsing command line arguments
	// -f: input fasta, -u: unit fasta, -s: score tsv, -v: variant tsv, -t: time tsv, -e: encodings txt
	int idx = 1;
	bool f_flag = false, u_flag = false, s_flag = false, v_flag = false, t_flag = false, e_flag = false;
	while(idx < argc){
		if(string(argv[idx]) == "-f" && idx + 1 < argc){
			read_file = argv[idx + 1];
			f_flag = true;
			idx += 2;
		}else if(string(argv[idx]) == "-u" && idx + 1 < argc){
			unit_file = argv[idx + 1];
			u_flag = true;
			idx += 2;
		}else if(string(argv[idx]) == "-s" && idx + 1 < argc){
			score_file = argv[idx + 1];
			s_flag = true;
			idx += 2;
		}else if(string(argv[idx]) == "-v" && idx + 1 < argc){
			variant_file = argv[idx + 1];
			v_flag = true;
			idx += 2;
		}else if(string(argv[idx]) == "-t" && idx + 1 < argc){
			time_file = argv[idx + 1];
			t_flag = true;
			idx += 2;
		}else if(string(argv[idx]) == "-e" && idx + 1 < argc){
			encodings_file = argv[idx + 1];
			e_flag = true;
			idx += 2;
		}else{
			cerr << "Usage: " << argv[0] << " -f <read_fasta> -u <unit_fasta> [-s <score_tsv>] [-v <variant_tsv>] [-t <time_txt>] [-e <encodings_txt>]" << endl;
			return 1;
		}
	}
	if(!f_flag || !u_flag){
		cerr << "Usage: " << argv[0] << " -f <read_fasta> -u <unit_fasta> [-s <score_tsv>] [-v <variant_tsv>] [-t <time_txt>] [-e <encodings_txt>]" << endl;
		return 1;
	}

	// Parsing input files (duplicated units are removed)
	vector<string> read_names;
	vector<vector<int>> reads;
	vector<vector<int>> units;
	read_fasta2(read_file, unit_file, read_names, reads, units);

	// Scores of mutations, insertions, deletions
	vector<vector<int>> eddc_units = {{0}, {1}, {2}, {3}};
	for(auto vec : units){eddc_units.push_back(vec);}
	Params params(1.0, 1.0, 0.5, eddc_units);

	// String decomposer
	vector<vector<int>> encodings;
	vector<vector<vector<int>>> decomposed_seqs;
	string_decompose(reads, read_names, units, encodings, decomposed_seqs, encodings_file, e_flag);

	// f score calculation limit length (setting minimum value of the limit)
	int f_par = 5
	int ulen_max = 0, ulen_min = INT_MAX;
	for(auto &unit : units){
		ulen_max = max(ulen_max, (int)unit.size());
		ulen_min = min(ulen_min, (int)unit.size());
	}
	f_par = max(f_par, ulen_max/ulen_min + 1);

	// Execute heddc_acc
	vector<vector<Score>> scores;
	vector<long long> measure_time;
	auto start = chrono::high_resolution_clock::now();
	heddc_all(encodings, units, params, scores, f_par, measure_time);
	auto end = chrono::high_resolution_clock::now();
	auto dur_msec = chrono::duration_cast<chrono::milliseconds>(end - start).count();
	auto dur_usec = chrono::duration_cast<chrono::microseconds>(end - start).count();

	// Results output
	if(s_flag) out_scores(scores, reads, score_file);
	else print_scores(scores, reads);
	if(v_flag) out_variants(scores, variant_file);
	if(t_flag) out_time(dur_msec, dur_usec, time_file, measure_time);

	return 0;
}