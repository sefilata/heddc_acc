#include "string_decomposer.h"


SDParams::SDParams() : match(1.0), mismatch(-1.0), indel(-1.0) {}
SDParams::SDParams(double m, double mis, double in) : match(m), mismatch(mis), indel(in) {}


SD::SD(int seq_len, const vector<vector<int>> &blocks, const SDParams &params){
	// resize
	dp.resize(blocks.size(), vector<vector<double>>());
	for(int b = 0; b < blocks.size(); b++){
		dp[b] = vector<vector<double>>(blocks[b].size()+1, vector<double>(seq_len+1, 0));
	}
	block_switch.resize(seq_len+1, -1);

	// initialize
	for(int j = 1; j < seq_len+1; j++){
		// for(int b = 0; b < blocks.size(); b++){dp[b][0][j] = j * params.indel;}
		for(int b = 0; b < blocks.size(); b++){dp[b][0][j] = -DBL_MAX;}
	}
}


void string_decomposer(
	const vector<int> &seq,
	const vector<vector<int>> &blocks,
	const SDParams &params,
	vector<int> &encoding,					// encoded by block ids
	vector<vector<int>> &decomposed_seq		// encoded by block sequences (with mutations)
){
	SD sd(seq.size(), blocks, params);

	// main DP
	for(int j = 1; j < seq.size()+1; j++){
		for(int b = 0; b < blocks.size(); b++){
			for(int i = 1; i < blocks[b].size()+1; i++){
				double match_score = (blocks[b][i-1] == seq[j-1]) ? params.match : params.mismatch;
				sd.dp[b][i][j] = max({
					sd.dp[b][i-1][j-1] + match_score,
					sd.dp[b][i][j-1] + params.indel,
					sd.dp[b][i-1][j] + params.indel,
					sd.dp[b][0][j-1] + (i-1)*params.indel + match_score		// パラメータによってはindelのみの場合とも比較が必要かも
				});
			}
		}
		// block-switching
		double last_max = -DBL_MAX;
		for(int b = 0; b < blocks.size(); b++){
			if(sd.dp[b][blocks[b].size()][j] > last_max){
				last_max = sd.dp[b][blocks[b].size()][j];
				sd.block_switch[j] = b;
			}
		}
		for(int b = 0; b < blocks.size(); b++){
			sd.dp[b][0][j] = last_max;
			// if(j < seq.size()){sd.dp[b][0][j+1] = last_max + params.indel;}
		}
	}

	// debug display
	// cout << "DP table: " << endl;
	// for(int b = 0; b < blocks.size(); b++){
	// 	cout << "(block " << b << ")" << endl;
	// 	for(int i = 0; i < blocks[b].size()+1; i++){
	// 		for(int j = 0; j < seq.size()+1; j++){
	// 			cout << setprecision(2) << sd.dp[b][i][j] << " ";
	// 		}cout << endl;
	// 	}
	// }
	// cout << "block switch: " << endl;
	// for(int j = 0; j < seq.size()+1; j++){cout << sd.block_switch[j] << " ";}
	// cout << endl;

	// back tracking
	vector<int> bt = {0, 0, (int)seq.size()};

	while(bt[2] > 0){
		if(bt[1] == 0){
			bt[0] = sd.block_switch[bt[2]];
			bt[1] = blocks[bt[0]].size();
			encoding.push_back(bt[0]);
			decomposed_seq.push_back({});
		}
		// cout << "back tracking: " << bt[0] << ", " << bt[1] << ", " << bt[2] << endl;

		double match_score = (blocks[bt[0]][bt[1]-1] == seq[bt[2]-1]) ? params.match : params.mismatch;
		double dp_score = sd.dp[bt[0]][bt[1]][bt[2]];
		if(abs(dp_score - (sd.dp[bt[0]][bt[1]-1][bt[2]-1] + match_score)) < DBL_ERR){
			decomposed_seq[decomposed_seq.size()-1].push_back(seq[bt[2]-1]);
			bt = {bt[0], bt[1]-1, bt[2]-1};
		}else if(abs(dp_score - (sd.dp[bt[0]][bt[1]][bt[2]-1] + params.indel)) < DBL_ERR){
			decomposed_seq[decomposed_seq.size()-1].push_back(seq[bt[2]-1]);
			bt = {bt[0], bt[1], bt[2]-1};
		}else if(abs(dp_score - (sd.dp[bt[0]][bt[1]-1][bt[2]] + params.indel)) < DBL_ERR){
			bt = {bt[0], bt[1]-1, bt[2]};
		}else if(abs(dp_score - (sd.dp[bt[0]][0][bt[2]-1] + (bt[1]-1)*params.indel + match_score)) < DBL_ERR){
			decomposed_seq[decomposed_seq.size()-1].push_back(seq[bt[2]-1]);
			bt = {bt[0], 0, bt[2]-1};
		}else{
			cerr << "ERROR(string_decomposer): No matching back-tracking path exists." << endl;
		}
	}

	reverse(encoding.begin(), encoding.end());
	reverse(decomposed_seq.begin(), decomposed_seq.end());
	for(auto &vec : decomposed_seq){reverse(vec.begin(), vec.end());}
}


// 以下，ユニット分解（あまりうまくいっていない）
#pragma region unit_decompose

// decompがunitのprefixに誤差1-similarity未満で一致するかどうかを判定し，一致するならばunits_to_addに追加する
void match_prefix(
	const vector<int> &unit,
	const vector<int> &decomp,
	vector<vector<int>> &units_to_add,
	const double similarity,
	const double mut_score,
	const double indel_score
){
	int n = unit.size();
	int m = decomp.size();

	vector<vector<double>> dp(n+1, vector<double>(m+1, 0.0));
	for(int i = 0; i <= n; i++){dp[i][0] = i * indel_score;}
	for(int i = 0; i <= m; i++){dp[0][i] = i * indel_score;}

	for(int i = 1; i <= n; i++){for(int j = 1; j <= m; j++){
		double mut_ij = (unit[i-1] == decomp[j-1]) ? 0.0 : mut_score;
		dp[i][j] = min({
			dp[i-1][j-1] + mut_ij,
			dp[i-1][j] + indel_score,
			dp[i][j-1] + indel_score
		});
	}}

	int min_idx = 0;
	for(int i = 1; i < n+1; i++){
		if(dp[i][m] <= dp[min_idx][m]){min_idx = i;}
	}
	// cout << "dp[min_idx][m]: " << dp[min_idx][m] << endl;
	// cout << (1.0-similarity)*(1.0*similarity)*min_idx*m << endl;
	if(dp[min_idx][m]*dp[min_idx][m] < (1.0-similarity)*(1.0*similarity)*min_idx*m){
		vector<int> rotate(n);
		for(int i = 0; i < n; i++){rotate[i] = unit[(i+min_idx)%n];}
		units_to_add.push_back(rotate);
	}

	// debug
	// cout << "[debug] match_prefix::dp" << endl;
	// for(int i = 0; i < n+1; i++){for(int j = 0; j < m+1; j++){
	// 	cout << dp[i][j] << "  ";
	// }cout << endl;}
}

// decompがunitのsuffixに誤差1-similarity未満で一致するかどうかを判定し，一致するならばunits_to_addに追加する
void match_suffix(
	const vector<int> &unit,
	const vector<int> &decomp,
	vector<vector<int>> &units_to_add,
	const double similarity,
	const double mut_score,
	const double indel_score
){
	vector<int> _unit(unit.size()), _decomp(decomp.size());
	vector<vector<int>> _units_to_add;
	for(int i = 0; i < unit.size(); i++){_unit[i] = unit[unit.size()-i-1];}
	for(int i = 0; i < decomp.size(); i++){_decomp[i] = decomp[decomp.size()-i-1];}

	match_prefix(_unit, _decomp, _units_to_add, similarity, mut_score, indel_score);

	for(auto vec : _units_to_add){
		reverse(vec.begin(), vec.end());
		units_to_add.push_back(vec);
	}
}

// 2つのvector<int>が同じかどうか判定する
bool is_same_vec(const vector<int> &vec1, const vector<int> &vec2){
	if(vec1.size() != vec2.size()){return false;}
	for(int i = 0; i < vec1.size(); i++){
		if(vec1[i] != vec2[i]){return false;}
	}
	return true;
}

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
){
	if(total_sub > max_diff){return false;}
	for(int i = 0; i < y_vec.size(); i++){
		if(x_vec[i+x_idx] != y_vec[i+y_idx]){total_sub++;}
		if(total_sub > max_diff){return false;}
	}
	return true;
}

bool unit_dfs(
	const vector<vector<int>> &blocks,
	const int idx,
	const int st,
	int &total_sub,
	const int max_diff,
	vector<vector<int>> &decomposition
){
	if(total_sub > max_diff){return false;}
	for(int i = 0; i < blocks.size(); i++){
		if(i == idx){continue;}
		if(st + blocks[i].size() > blocks[idx].size()){continue;}
		int total_sub_cp = total_sub;
		if(sub_only_dp(blocks[idx], blocks[i], st, 0, total_sub_cp, max_diff)){
			decomposition.push_back({st, st+(int)blocks[i].size()});
			if(blocks[idx].size() == st+blocks[i].size() && total_sub <= max_diff){return true;}
			if(unit_dfs(blocks, idx, st+blocks[i].size(), total_sub_cp, max_diff, decomposition)){
				return true;
			}
			decomposition.pop_back();
		}
	}
	return false;
}

// ユニット分解
void unit_dfs_init(
	vector<vector<int>> &blocks,
	const double similarity,
	vector<vector<vector<int>>> &decompositions
){
	// vector<vector<vector<int>>> decompositions;		// [unit to decompose], [subunit 1], [subunit 2], ...
	for(int i = 0; i < blocks.size(); i++){
		// cout << "max_diff: " << endl;
		// cout << i << ": " << (int)((1-similarity) * blocks[i].size()) << endl;
		vector<vector<int>> decomposition = {{i}};
		int total_sub = 0;
		if(unit_dfs(blocks, i, 0, total_sub, (int)((1-similarity) * blocks[i].size()), decomposition)){
			// 分解すべきユニットに近い別のユニットがある場合，正しく分解されない場合がある。
			// <<要修正>>
			if(decomposition.size() <= 2){continue;}
			decompositions.push_back(decomposition);
		}
	}

	// debug
	// cout << "decompositions: " << endl;
	// for(auto dec : decompositions){
	// 	for(auto i : dec){
	// 		cout << "[";
	// 		for(int j : i){cout << j << ' ';}
	// 		cout << "]";
	// 	}cout << endl;
	// }
}

// vectorのスライスを取得
vector<int> vec_copy(const vector<int> &vec, int st, int ed){
	vector<int> ret(ed-st);
	for(int i = 0; i < ed-st; i++){ret[i] = vec[st+i];}
	return ret;
}

// blocks, encodings, decomposed_seqsを更新する
// 新しいユニットは既出でないかどうかだけ確認して追加する。古いユニットは削除しない。
void update_blocks(
	vector<vector<int>> &blocks,
	vector<vector<int>> &encodings,
	vector<vector<vector<int>>> &decomposed_seqs,
	const vector<vector<vector<int>>> &decompositions
){
	// ユニット置換の対応表
	map<int, vector<int>> unit_map;

	// 対応表の作成とblocksの更新
	for(auto dec : decompositions){
		int idx = dec[0][0];		// 置換元のユニット番号
		vector<int> new_units(0);	// 置換先のユニット番号列
		for(int i = 1; i < dec.size(); i++){
			int st = dec[i][0];
			int ed = dec[i][1];
			int in_blocks = -1;		// すでにblocksに含まれているかどうか？
			for(int block = 0; block < blocks.size(); block++){
				if(blocks[block].size() != ed - st){continue;}
				bool to_continue = false;
				for(int j = 0; j < blocks[block].size(); j++){
					if(blocks[block][j] != blocks[idx][st+j]){
						to_continue = true;
						break;
					}
				}
				if(to_continue){continue;}
				// cout << "dec in a block: " << idx << ", st=" << st << ", ed=" << ed << ", " << block << endl;
				in_blocks = block;
			}
			if(in_blocks < 0){
				blocks.push_back(vec_copy(blocks[idx], st, ed));
				new_units.push_back(blocks.size()-1);
			}else{
				new_units.push_back(in_blocks);
			}
		}
		unit_map[idx] = new_units;
	}

	// debug
	// cout << "unit_map: " << endl;
	// for(const auto& [key, vec] : unit_map){
	// 	cout << key << ": ";
	// 	for(auto i : vec){cout << i << ' ';}
	// 	cout << endl;
	// }

	// encodingsの更新
	for(int i = 0; i < encodings.size(); i++){
		vector<int> new_vec(0);
		for(int j : encodings[i]){
			if(unit_map.find(j) == unit_map.end()){
				new_vec.push_back(j);
			}else{
				for(int k : unit_map[j]){new_vec.push_back(k);}
			}
		}
		encodings[i] = new_vec;
	}

	// decomposed_seqsの更新
	// （省略）
}

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
	const double similarity,		// unitの分解時の類似度の下限値
	const int max_diff		// 連続する何塩基の誤差まで許容するか？
){
	// encodingでのblockの境を特定する
	// 境の前側のブロックを{seq番号, ブロックの順番} (両方0-indexed)で記録する
	vector<pair<int,int>> borders;
	for(int i = 0; i < encodings.size(); i++){
		int prev = encodings[i][0];
		for(int j = 1; j < encodings[i].size(); j++){
			if(prev == encodings[i][j]){continue;}
			borders.emplace_back(i, j-1);
			prev = encodings[i][j];
		}
	}

	// debug
	// cout << "encodings: " << endl;
	// for(auto vec : encodings){
	// 	for(int i : vec){cout << i << ' ';}
	// 	cout << endl;
	// }
	// cout << "borders: " << endl;
	// for(auto p : borders){
	// 	cout << "[" << p.first << ", " << p.second << "] :  ";
	// 	for(int i : decomposed_seqs[p.first][p.second]){cout << i;}
	// 	cout << ", ";
	// 	for(int i : decomposed_seqs[p.first][p.second+1]){cout << i;}
	// 	cout << endl;
	// }

	// bordersの各decomposed_seqがprefixおよびsuffixと一致するかどうか確認し，一致すれば追加ユニット候補とする
	// 最初および最後のユニットは考慮しない
	vector<vector<int>> units_to_add;
	vector<int> units_to_del;
	for(auto [f, s] : borders){
		// 境界の前方
		int units_num = units_to_add.size();
		if(decomposed_seqs[f][s].size() <= blocks[encodings[f][s]].size() && s > 0){
			match_prefix(blocks[encodings[f][s]], decomposed_seqs[f][s], units_to_add, similarity);
		}
		if(units_num != units_to_add.size()){
			units_num = units_to_add.size();
			units_to_del.push_back(encodings[f][s]);
		}

		// 境界の後方
		if(decomposed_seqs[f][s+1].size() > blocks[encodings[f][s+1]].size() && s+1 < decomposed_seqs[f].size()-1){
			match_suffix(blocks[encodings[f][s+1]], decomposed_seqs[f][s+1], units_to_add, similarity);
		}
		if(units_num != units_to_add.size()){
			units_to_del.push_back(encodings[f][s+1]);
		}
	}

	// 回転前のunitを削除
	sort(units_to_del.rbegin(), units_to_del.rend());
	auto last = std::unique(units_to_del.begin(), units_to_del.end());
	units_to_del.erase(last, units_to_del.end());
	// cout << "units_to_del: " << endl;
	// for(int i : units_to_del){cout << i << ' ';}
	// cout << endl;
	for(int i : units_to_del){blocks.erase(blocks.begin() + i);}

	// 重複を除いてから再度string decomposer
	sort(units_to_add.begin(), units_to_add.end());
	units_to_add.erase(unique(units_to_add.begin(), units_to_add.end()), units_to_add.end());
	int unit_num = blocks.size();
	for(auto vec : units_to_add){
		bool in_blocks = false;
		for(int i = 0; i < unit_num; i++){
			if(is_same_vec(blocks[i], vec)){
				in_blocks = true;
				break;
			}
		}
		if(in_blocks){continue;}
		blocks.push_back(vec);
	}
	for(int i = 0; i < reads.size(); i++){
		encodings[i].clear();
		decomposed_seqs[i].clear();
		string_decomposer(reads[i], blocks, params, encodings[i], decomposed_seqs[i]);
	}

	// debug
	// cout << "blocks: " << endl;
	// for(auto vec : blocks){
	// 	for(int i : vec){cout << i << ' ';}
	// 	cout << endl;
	// }
	// cout << "decomposed_seqs: " << endl;
	// for(auto vec : decomposed_seqs){
	// 	for(auto v : vec){
	// 		for(int i : v){cout << i;}
	// 		cout << ' ';
	// 	}
	// 	cout << endl;
	// }

	// ユニットの分解
	vector<vector<vector<int>>> decompositions;
	unit_dfs_init(blocks, similarity, decompositions);
	update_blocks(blocks, encodings, decomposed_seqs, decompositions);
}

// unit variant を指定した個数 (-1の場合は全て) ユニットとして追加して再度 string decomposer かける
void unit_variant(
	const vector<vector<int>> &reads,
	vector<vector<int>> &blocks,
	const SDParams &params,
	vector<vector<int>> &encodings,
	vector<vector<vector<int>>> &decomposed_seqs,
	int num_variants		// 何個のvariantsをとるか？
){
	if(num_variants != -1){
		cerr << "Warning: num_variants is not yet implemented." << endl;
	}

	set<vector<int>> variant_set;
	for(const auto &enc_list : decomposed_seqs){
		for(const auto &encoding : enc_list){
			variant_set.insert(encoding);
		}
	}
	// blocksの更新 (上書き)
	blocks.assign(variant_set.begin(), variant_set.end());

	// 再度 string decomposer
	for(int i = 0; i < reads.size(); i++){
		encodings[i].clear();
		decomposed_seqs[i].clear();
		string_decomposer(reads[i], blocks, params, encodings[i], decomposed_seqs[i]);
	}
}

#pragma endregion
