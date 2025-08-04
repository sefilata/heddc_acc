#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <math.h>
#include <chrono>

#include "eddc.h"
#include "fasta.h"
#include <omp.h>

using namespace std;


// csv形式で標準出力
void outResult(const vector<vector<double>> &trED){
	int read_num = trED.size();
	for(int i = 0; i < read_num; i++){
		for(int j = 0; j < read_num; j++){
			cout << trED[i][j];
			if(j != read_num - 1){cout << ", ";}
		}
		cout << endl;
	}
}

int main(int argc, char* argv[]){
	// TR配列
	// vector<int> s = {0,1,0,1,0,1};
	// vector<int> t = {0,1,0,1,0,1,0,1};
	// vector<int> s = {0,1,1,0,1,1,0,1,1,3,2,3,2,3,2};
	// vector<int> t = {0,1,1,0,1,1,3,2,3,2,3,2};

	// unit: 最初4つはbase(0, 1, 2, 3)で固定
	// vector<vector<int>> units = {{0}, {1}, {2}, {3}, {0,1,1}, {3,2}, {0,1}};
	// vector<vector<int>> units = {{0}, {1}, {2}, {3}, {0,1}};

	string read_file = argv[1];
	string unit_file = argv[2];
	vector<string> read_names;
	vector<vector<int>> reads;
	vector<vector<int>> units = {{0}, {1}, {2}, {3}};

	read_fasta2(read_file, unit_file, read_names, reads, units);

	// 確認
	// cout << reads.size() << endl;
	// for(int i = 0; i < read_names.size(); i++){
	// 	cout << read_names[i] << ": " << endl;
	// 	for(int b : reads[i]){cout << b;}
	// 	cout << endl;
	// }
	// cout << endl;
	// for(int i = 0; i < units.size(); i++){
	// 	for(int b : units[i]){cout << b;}
	// 	cout << endl;
	// }
	// cout << endl;

	// mut, indel, dupのスコア
	Params params(1.0, 1.0, 0.5, units);

	// // unit間編集距離の確認
	// for(int i = 0; i < units.size(); i++){
	// 	for(int j = 0; j < units.size(); j++){
	// 		cout << params.get_unit_to_unit(i, j).get_score() << "  ";
	// 	}
	// 	cout << endl;
	// }

	// EDDCの計算
	vector<vector<double>> scores(reads.size(), vector<double>(reads.size(), 0.0));
	// Score val = calc_eddc(s, t, units, params);
	// Score val = calc_eddc(reads[0], reads[1], units, params);
	int n = reads.size();
	omp_set_num_threads(30);
	#pragma omp parallel for schedule(dynamic)
	for(int s = 0; s < n*n; s++){
		int i = s / n;
		int j = s % n;
		if(i >= j){continue;}
		// #pragma omp critical
		// cerr << i << ", " << j << endl;

		Score val = calc_eddc(reads[i], reads[j], units, params);
		double normalized = val.get_score() / sqrt(reads[i].size() * reads[j].size());

		#pragma omp critical
		{
			scores[i][j] = normalized;
			scores[j][i] = normalized;
			cerr << "calculated the distance between " << i << " and " << j << endl;
		}
	}
	cerr << "finished calculating all scores" << endl;

	// 結果の表示
	// cout << "score: " << val.get_score() << endl;
	// cout << "mut: " << val.get_mut() << ", indel: " << val.get_indel() << ", dup: [";
	// val.print_dup();
	// cout << "]" << endl;

	// 結果の書き込み
	outResult(scores);
}