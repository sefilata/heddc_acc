#pragma once
#include <iostream>
#include <vector>
#include <algorithm>
#include <cfloat>
#include <iomanip>
#include <map>
#include <set>

#define DBL_ERR 1e-10

using namespace std;

class SDParams{
public:
	const double match;
	const double mismatch;
	const double indel;

	SDParams();
	SDParams(double m, double mis, double in);
};

class SD{
public:
	vector<vector<vector<double>>> dp;
	vector<int> block_switch;		// back trackingで使用

	SD(int seq_len, const vector<vector<int>> &blocks, const SDParams &params);
};

void string_decomposer(
	const vector<int> &seq,
	const vector<vector<int>> &blocks,
	const SDParams &params,
	vector<int> &encoding,					// encoded by block ids
	vector<vector<int>> &decomposed_seq		// encoded by block sequences (with mutations)
);


// decompがunitのprefixに誤差1-similarity未満で一致するかどうかを判定し，一致するならばunits_to_addに追加する
void match_prefix(
	const vector<int> &unit,
	const vector<int> &decomp,
	vector<vector<int>> &units_to_add,
	const double similarity,
	const double mut_score = 1.0,
	const double indel_score = 1.0
);


// decompがunitのsuffixに誤差1-similarity未満で一致するかどうかを判定し，一致するならばunits_to_addに追加する
void match_suffix(
	const vector<int> &unit,
	const vector<int> &decomp,
	vector<vector<int>> &units_to_add,
	const double similarity,
	const double mut_score = 1.0,
	const double indel_score = 1.0
);

// 2つのvector<int>が同じかどうか判定する
bool is_same_vec(const vector<int> &vec1, const vector<int> &vec2);


// x_idx, y_idxからDPを開始する
// bool diagonal_dp(
// 	const vector<int> &x_vec,
// 	const vector<int> &y_vec,
// 	int x_idx,
// 	int y_idx,
// 	vector<int> &dp_ini,		// 初期値 (最後に更新する)
// 	const double similarity,	// 類似度の下限値
// 	const int max_diff			// 連続する何塩基の誤差まで許容するか？
// 	// 連続する何塩基というのではこの方法は使えない。全体で何塩基というなら使えるが、、
// ){
// 	int width = max_diff*2 + 1;
// 	vector<vector<int>> dp(2, vector<int>(width, 0));
// }


// indelを考慮しない場合
bool sub_only_dp(
	const vector<int> &x_vec,
	const vector<int> &y_vec,
	int x_idx,
	int y_idx,
	int &total_sub,
	const int max_diff
);


bool unit_dfs(
	const vector<vector<int>> &blocks,
	const int idx,
	const int st,
	int &total_sub,
	const int max_diff,
	vector<vector<int>> &decomposition
);

// ユニット分解
void unit_dfs_init(
	vector<vector<int>> &blocks,
	const double similarity,
	vector<vector<vector<int>>> &decompositions
);

// vectorのスライスを取得
vector<int> vec_copy(const vector<int> &vec, int st, int ed);

// blocks, encodings, decomposed_seqsを更新する
// 新しいユニットは既出でないかどうかだけ確認して追加する。古いユニットは削除しない。
void update_blocks(
	vector<vector<int>> &blocks,
	vector<vector<int>> &encodings,
	vector<vector<vector<int>>> &decomposed_seqs,
	const vector<vector<vector<int>>> &decompositions
);

// ユニットが切り替わるときの前後のdecomposed_seqを見て，それらが理想的なユニットのprefixやsuffixになっていないか調べる。
// そうなっていた場合，その点で切って回転させたユニットを追加して再度string decomposerにかける。
// その後各ユニットに対して他のより短いユニットの組み合わせで書けないかどうか判定する。
// ← 貪欲に組み立てていく？
void unit_decompose(
	const vector<vector<int>> &reads,
	vector<vector<int>> &blocks,
	const SDParams &params,
	vector<vector<int>> &encodings,					// encoded by block ids
	vector<vector<vector<int>>> &decomposed_seqs,	// encoded by block sequences (with mutations)
	const double similarity = 0.8,		// unitの分解時の類似度の下限値
	const int max_diff = 1		// 連続する何塩基の誤差まで許容するか？
);


// unit variant を指定した個数 (-1の場合は全て) ユニットとして追加して再度 string decomposer かける
void unit_variant(
	const vector<vector<int>> &reads,
	vector<vector<int>> &blocks,
	const SDParams &params,
	vector<vector<int>> &encodings,
	vector<vector<vector<int>>> &decomposed_seqs,
	int num_variants = -1		// 何個のvariantsをとるか？
);
