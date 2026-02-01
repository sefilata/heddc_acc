#pragma once
#include <iostream>
#include <vector>
#include <chrono>
#include <map>
#include <set>
#include <algorithm>
#include <iomanip>
#include <cmath>
#include <cfloat>
#include <tuple>
#include <compare>
#include <chrono>
#include <climits>
#include <fstream>
#include <string>
#include <omp.h>
#include <numeric>
#include <sstream>

#include "../string_decomposer/string_decomposer.h"
#include "../eddc_original/fasta.h"

using namespace std;

const bool MEASURE_TIME = false;
const bool MEASURE_MEMORY_LINUX = true;
const bool DEBUG_OUT = false;
const double PARAM_0_VALUE = 0.5; 	// 0: 0.5*unit_length
const double PARAM_1_VALUE = 0.6; 	// 1: all

class Params;

// Score class
class Score{
private:
	double score;
	int mut, indel;
	map<int, int> dup;

public:
	Score();
	Score(double s);

	// Adding mutations
	void set_indel(int len, const Params &par);
	// void set_mut(const Params &par);
	// void set_dup(int unit_id, const Params &par);
	void set_zero();

	// Getting variables
	double get_score() const;
	int get_mut() const;
	int get_indel() const;
	void print_dup(ostream& os, const vector<vector<int>> &units) const;

	// Adding mutations (used only in edit_distance() and Params constructor)
	void set_indel_2(int len, double par);
	void set_mut_2(double par);
	void set_dup_2(int unit_id, double par);

	// Addition of Score
	Score operator+(const Score &other) const;
	Score& operator+=(const Score &other);

	// Order of Score
	auto operator<=>(const Score &other) const;
	bool operator==(const Score &other) const;
};

// Conventional edit distance
Score edit_distance(const vector<int> &s, const vector<int> &t, double mut, double indel, int st = 0, int ed = -1);

// Params class functions
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

	// Adding mutations
	void set_indel(int len, const Params &par);
	void set_mut(const Params &par);
	void set_dup(int unit_id, const Params &par);
};

// DEBUG
void print_vec2_score(const vector<vector<Score>> &vec2, const vector<vector<int>> &units);

// Measure memory (Linux only)
long long get_hwm_kb();


/* ---------------- Naive EDDC Starts Here ---------------- */

// Stage1 (Almost the same as in the original paper)
// The computations for s and t are the same
// units: [{bases}. {units}, eps]
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

/* ---------------- Naive EDDC Ends Here ---------------- */


// heddc(s[i,j), x)
void calc_f(
	// reads は encoded reads のこと（あとで修正）
	const vector<vector<int>> &reads,
	const vector<vector<int>> &units,
	const Params &params,
	vector<vector<vector<vector<Score>>>> &f_scores,
	const vector<int> &acc_par1s,
	vector<long long> &measure_time,
	vector<long long> &measure_mem
);

// Compute hEDDC given two TRs and f_scores
Score calc_heddc(
	// u and v are unit-encoded sequences
	const vector<int> &u,	// seq1
	int u_idx,
	const vector<int> &v,	// seq2
	int v_idx,
	const vector<vector<int>> &units,
	const Params &params,
	const vector<vector<vector<vector<Score>>>> &f_scores,
	vector<int> acc_par1s, 	// length limits when computing f_score
	int acc_par2 	// length limits when computing main dp (currently not used)
);

// Compute hEDDC for all TR pairs
void heddc_all(
	const vector<vector<int>> &encodings,
	vector<vector<int>> &units,
	const Params &params,
	vector<vector<Score>> &scores,
	int par,
	vector<long long> &measure_time,
	vector<long long> &measure_mem,
	int parallel_num_1,
	int parallel_num_2,
	vector<int> &tr_nums
);
