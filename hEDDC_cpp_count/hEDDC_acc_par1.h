#pragma once
#include <iostream>
#include <vector>
#include <chrono>
#include <map>
#include <set>
#include <algorithm>
#include <iomanip>
#include <math.h>
#include <cfloat>
#include <tuple>
#include <compare>
#include <chrono>
#include <climits>
#include <fstream>
#include <string>

#include "../string_decomposer/string_decomposer.h"
#include "../eddc_original/fasta.h"

using namespace std;

const bool MEASURE_TIME = false;

class Params;

// スコア構造体
class Score{
private:
	double score;
	int mut, indel;
	map<int, int> dup;

public:
	Score();
	Score(double s);

	// 各変数の追加
	void set_indel(int len, const Params &par);
	// void set_mut(const Params &par);
	// void set_dup(int unit_id, const Params &par);
	void set_zero();

	// 変数の取得
	double get_score() const;
	int get_mut() const;
	int get_indel() const;
	void print_dup(ofstream &ofs, const vector<vector<int>> &units) const;

	// edit_distance()とParamsのコンストラクタでのみ使用
	void set_indel_2(int len, double par);
	void set_mut_2(double par);
	void set_dup_2(int unit_id, double par);

	// Score同士の加算
	Score operator+(const Score &other) const;

	// Score同士の加算代入
	Score& operator+=(const Score &other);

	// Score同士の順序
	auto operator<=>(const Score &other) const;
};

// 従来の編集距離
Score edit_distance(const vector<int> &s, const vector<int> &t, double mut, double indel, int st = 0, int ed = -1);

// パラメータ構造体
// get_unit_to_unit() はこのままで良い？？
class Params{
public:
	double mut, indel;
	vector<vector<Score>> unit_to_unit;
	vector<Score> dup_scores;
	vector<Score> indel_scores;

	Params(double m, double i, int dup_pattern, const vector<vector<int>> &units);

	Score get_unit_to_unit(int i, int j) const;
	Score get_dup(int a) const;
	Score get_indel(int a) const;

	double mut_val() const;
	double indel_val() const;
	// double dup_val() const;

	// 変異の追加
	void set_indel(int len, const Params &par);
	void set_mut(const Params &par);
	void set_dup(int unit_id, const Params &par);
};

// unitsにbase (A,T,G,C)とepsilonを追加する
void add_base_eps(vector<vector<int>> &units);

// debug
void print_vec2(const vector<vector<int>> &vec2);
void print_vec2_score(const vector<vector<Score>> &vec2);
void print_vec3_score(vector<vector<vector<Score>>> vec);


/* ---------------- Naive EDDC ココカラ ---------------- */

// 以下 "../eddc_original/eddc.h" からコピーして少し改変したもの
// Stage1の計算(だいたい論文通り)
// sとtは方向が逆なだけなので計算方法は同じ
void string_to_unit(
	const vector<int> &s,
	const vector<vector<int>> &units,
	vector<vector<Score>> &S_eps,
	vector<vector<vector<Score>>> &S,
	vector<vector<vector<Score>>> &S2,
	const Params &params
);

Score calc_eddc(
	const vector<int> &s,
	const vector<int> &t,
	const vector<vector<int>> &units,
	const Params &params
);

/* ---------------- Naive EDDC ココマデ ---------------- */


// eddc_unitsの作成（仮）（あとで共通化する）
void shift_units(const vector<vector<int>> &units, vector<vector<int>> &eddc_units);

void calc_f(
	// reads は encoded reads のこと（あとで修正）
	const vector<vector<int>> &reads,
	const vector<vector<int>> &units,
	const Params &params,
	vector<vector<vector<vector<Score>>>> &f_scores,
	vector<int> acc_par1s,
	vector<long long> &measure_time
);

// 2本のTRとf_scoresが与えられたときにEDDCを計算する
Score calc_heddc(
	// u, v は encoded read
	const vector<int> &u,	// リード1
	int u_idx,
	const vector<int> &v,	// リード2
	int v_idx,
	const vector<vector<int>> &units,
	const Params &params,
	const vector<vector<vector<vector<Score>>>> &f_scores,
	vector<int> acc_par1s,		// f_score計算時の最大長さ
	int acc_par2 	// main dpでの最大長さ
);

// 全てのTRの組に対してEDDCを計算する
void heddc_all(
	const vector<vector<int>> &encodings,
	vector<vector<int>> &units,
	const Params &params,
	vector<vector<Score>> &scores,
	int par,
	vector<long long> &measure_time
);
