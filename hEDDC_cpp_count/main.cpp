#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <chrono>
#include <cstdlib>
#include <climits>
#include <iomanip>

#include "../eddc_original/fasta.h"
#include "heddc_acc_parallel.h"
#include "../string_decomposer/string_decomposer.h"

using namespace std;


int PARAMETER_PATTERN = 1; 	// Pattern of parameters of dup/cont ->  0: 0.5*unit_length, 1: all 0.6 (default)
int ACC_PATTERN = 1; 	// hEDDC or heddc_acc ->  0: hEDDC, 1: heddc_acc (default)
const int SCORE_PRECISION = 8; 	// Number of decimal places for score output
const int MIN_F_PAR = 5; 	// Minimum length limit for heddc_acc
const double MUTATION_COST = 1.0; 	// Mutation cost
const double INDEL_COST = 1.0; 	// Insertion/deletion cost


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

	// Not to include bases in string decomposer
	// Fix later
	for(auto &enc : encodings){
		for(auto &b : enc){b += 4;}
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
			ofs << setprecision(SCORE_PRECISION) << scores[i][j].get_score() / sqrt((double)reads[i].size() * (double)reads[j].size());;
			if(j != read_num-1) ofs << "\t";
		}
		ofs << endl;
	}
	ofs.close();
}

// Score output in stdout
void print_scores(
	const vector<vector<Score>> &scores,
	const vector<vector<int>> &reads
){
	int read_num = reads.size();

	for(int i = 0; i < read_num; i++){
		for(int j = 0; j < read_num; j++){
			cout << setprecision(SCORE_PRECISION) << scores[i][j].get_score() / sqrt((double)reads[i].size() * (double)reads[j].size());;
			if(j != read_num-1) cout << "\t";
		}
		cout << endl;
	}
}

// Output variations (mutations, indels, dupcons)
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

	for(int i = 0; i < (int)scores.size(); i++){
		for(int j = 0; j < (int)scores[i].size(); j++){
			ofs << "{mut:" << scores[i][j].get_mut();
			ofs << ", indel:" << scores[i][j].get_indel();
			ofs << ", dup:";
			scores[i][j].print_dup(ofs, units);
			ofs << "}";
			if(j != (int)scores[i].size()-1) ofs << "\t";
		}
		ofs << endl;
	}
	ofs.close();
}

// Output execution time
void out_time(
	long long dur_msec,
	long long dur_usec,
	string out_file,
	const vector<long long> &measure_time,
	const vector<long long> &measure_mem,
	int parallel_num_1,
	int parallel_num_2,
	const vector<int> &tr_nums
){
	ofstream ofs(out_file);
	if(!ofs.is_open()){
		cerr << "Error opening time output file: " << out_file << endl;
		return;
	}
	ofs << tr_nums[0] << " reads before compression" << endl;
	ofs << tr_nums[1] << " reads after compression" << endl;
	ofs << parallel_num_1 << " parallel (other than main dp)" << endl;
	ofs << parallel_num_2 << " parallel (main dp)" << endl;
	ofs << ACC_PATTERN << " ACC_PATTERN (0:hEDDC, 1:heddc_acc)" << endl;
	ofs << PARAMETER_PATTERN << " PARAMETER_PATTERN (0:0.5*|unit|, 1:all 0.6)" << endl;
	ofs << dur_msec << " msec" << endl;
	ofs << dur_usec << " usec" << endl;
	ofs << measure_time[0] << " msec (c1,c2,valid_rules)" << endl;
	ofs << measure_time[1] << " msec (f_scores)" << endl;
	ofs << measure_time[2] << " msec (main dp)" << endl;
	if(MEASURE_MEMORY_LINUX){
		ofs << measure_mem[0] << " kB MAX (after c1,c2, valid_rules)" << endl;
		ofs << measure_mem[1] << " kB MAX (after f_scores)" << endl;
		ofs << measure_mem[2] << " kB MAX (after main dp)" << endl;
	}
	ofs.close();
}

int main(int argc, char* argv[]){
	string read_file; 		// TR fasta (in)
	string unit_file; 		// Unit fasta (in)
	string score_file; 		// Score file (out)
	string variant_file; 	// Variations file (out)
	string time_file; 		// Execution time file (out)
	string encodings_file; 	// Etring decomposition results file (out)
	int parallel_num_1 = 1;	// The number of threads (other than main dp)
	int parallel_num_2 = 1;	// The number of threads (main dp)

	// Parsing command line arguments
	// -f: input fasta, -u: unit fasta, -s: score tsv, -v: variant tsv, -t: time tsv, -e: encodings txt, -p1: threads-1, -p2:threads-2 -p: parameter pattern, -a; acc pattern
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
		}else if(string(argv[idx]) == "-p1" && idx + 1 < argc){
			parallel_num_1 = atoi(argv[idx + 1]);
			if(parallel_num_1 < 1){
				cerr << "The number of threads must be at least 1." << endl;
				return 1;
			}
			idx += 2;
		}else if(string(argv[idx]) == "-p2" && idx + 1 < argc){
			parallel_num_2 = atoi(argv[idx + 1]);
			if(parallel_num_2 < 1){
				cerr << "The number of threads must be at least 1." << endl;
				return 1;
			}
			idx += 2;
		}else if(string(argv[idx]) == "-p" && idx + 1 < argc){
			PARAMETER_PATTERN = atoi(argv[idx + 1]);
			if(PARAMETER_PATTERN != 0 && PARAMETER_PATTERN != 1){
				cerr << "Parameter pattern must be 0 (0.5*|unit|) or 1 (all 0.6, default)." << endl;
				return 1;
			}
			idx += 2;
		}else if(string(argv[idx]) == "-a" && idx + 1 < argc){
			ACC_PATTERN = atoi(argv[idx + 1]);
			if(ACC_PATTERN != 0 && ACC_PATTERN != 1){
				cerr << "Acc pattern must be 0 (hEDDC) or 1 (heddc_acc, default)." << endl;
				return 1;
			}
			idx += 2;
		}else{
			cerr << "Usage: " << argv[0] << " -f <read_fasta> -u <unit_fasta> [-s <score_tsv>] [-v <variant_tsv>] [-t <time_txt>] [-e <encodings_txt>] [-p1 <parallel_num_1>] [-p2 <parallel_num_2>] [-p <parameter_pattern>] [-a <acc_pattern>]" << endl;
			return 1;
		}
	}
	if(!f_flag || !u_flag){
		cerr << "Usage: " << argv[0] << " -f <read_fasta> -u <unit_fasta> [-s <score_tsv>] [-v <variant_tsv>] [-t <time_txt>] [-e <encodings_txt>] [-p1 <parallel num (f-scores)>] [-p2 <parallel num (main-dp)>] [-p <parameter_pattern>] [-a <acc_pattern>]" << endl;
		return 1;
	}


	// Parsing input files (duplicated units are removed)
	vector<string> read_names;
	vector<vector<int>> reads;
	vector<vector<int>> units_only; 	// units (without bases and epsilon)
	read_fasta2(read_file, unit_file, read_names, reads, units_only);

	// Scores of mutations, insertions, deletions
	vector<vector<int>> units = {{0}, {1}, {2}, {3}}; 	// single bases (A,T,G,C)
	for(auto u : units_only){units.push_back(u);}
	units.push_back({}); 	// empty string (epsilon)
	Params params(MUTATION_COST, INDEL_COST, PARAMETER_PATTERN, units);

	// String decomposer
	vector<vector<int>> encodings;
	vector<vector<vector<int>>> decomposed_seqs;
	string_decompose(reads, read_names, units_only, encodings, decomposed_seqs, encodings_file, e_flag);

	// Computation range (length limit) for heddc_acc
	// (iii) Case: Units with substantially different lengths.
	int f_par = MIN_F_PAR; 	// Minimum length limit
	if(ACC_PATTERN == 0){
		f_par = INT_MAX; 	// hEDDC
	}else{
		int ulen_max = 0, ulen_min = INT_MAX;
		for(int i = 4; i < (int)units.size()-1; i++){ 	// excluding single bases and epsilon
			ulen_max = max(ulen_max, (int)units[i].size());
			ulen_min = min(ulen_min, (int)units[i].size());
		}
		f_par = max(f_par, ulen_max/ulen_min + 1);
	}

	// EDDC
	vector<vector<Score>> scores;
	vector<long long> measure_time(3, -1LL);
	vector<long long> measure_mem(3, -1LL);
	vector<int> tr_nums(2, -1); 	// trnums[0]: # of seqs before compression, trnums[1]: # of seqs after compression
	auto start = chrono::high_resolution_clock::now();
	heddc_all(encodings, units, params, scores, f_par, measure_time, measure_mem, parallel_num_1, parallel_num_2, tr_nums);
	auto end = chrono::high_resolution_clock::now();
	auto dur_msec = chrono::duration_cast<chrono::milliseconds>(end - start).count();
	auto dur_usec = chrono::duration_cast<chrono::microseconds>(end - start).count();

	// Results output
	if(s_flag) out_scores(scores, reads, score_file);
	else print_scores(scores, reads);
	if(v_flag) out_variants(scores, variant_file, units);
	if(t_flag) if(t_flag) out_time(dur_msec, dur_usec, time_file, measure_time, measure_mem, parallel_num_1, parallel_num_2, tr_nums);

	if(DEBUG_OUT){
		auto now = chrono::system_clock::now();
		time_t t = chrono::system_clock::to_time_t(now);
		cout << "Finished all at " << ctime(&t);
	}

	return 0;
}