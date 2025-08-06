#include "hEDDC_acc_par1.h"


// Constructor of Score
Score::Score() : score(0), mut(0), indel(0), dup() {}
Score::Score(double s) : score(s), mut(0), indel(0), dup() {}
// Setting zero for Score
void Score::set_zero(){score = 0;}
// Getting variables in Score
double Score::get_score() const {return score;}
int Score::get_mut() const {return mut;}
int Score::get_indel() const {return indel;}
void Score::print_dup(ofstream &ofs, const vector<vector<int>> &units) const {
	for(const auto &dup_ : dup){
		if(dup_.second > 0){
			if(dup_.first < 4){
				cerr << "Warning: print_dup() called with unit_id < 4, which is not expected." << endl;
				continue;
			}
			ofs << "(";
			string unit_str;
			digits_to_string(units[dup_.first - 4], unit_str);
			ofs << unit_str;
			ofs << ", " << dup_.second << ")";
		}
	}
}
// Used only in edit_distance() and Constructor of Params
void Score::set_indel_2(int len, double par){score += par * len; indel += len;}
void Score::set_mut_2(double par){score += par; mut++;}
void Score::set_dup_2(int unit_id, double par){score += par; dup.emplace(unit_id, 1);}
// Addition of Scores
Score Score::operator+(const Score &other) const{
	if(score < 0 || other.score < 0){cerr << "Error: operator + for minus scores" << endl;}
	Score ret(score);
	ret.score += other.score;
	ret.mut = mut + other.mut;
	ret.indel = indel + other.indel;
	ret.dup.insert(dup.begin(), dup.end());
	for(auto it = other.dup.begin(); it != other.dup.end(); it++){
		if(ret.dup.find(it->first) == ret.dup.end()){
			ret.dup.insert(*it);
		}else{
			ret.dup[it->first] += it->second;
		}
	}
	return ret;
}
// Addition assignment of Scores
Score& Score::operator+=(const Score &other){
	if(score < 0 || other.score < 0){cerr << "Error: operator += for minus scores" << endl;}
	score += other.score;
	mut += other.mut;
	indel += other.indel;
	for(auto it = other.dup.begin(); it != other.dup.end(); it++){
		if(dup.find(it->first) == dup.end()){
			dup.insert(*it);
		}else{
			dup[it->first] += it->second;
		}
	}
	return *this;
}
// Comparison of Scores
auto Score::operator<=>(const Score &other) const{
	return score <=> other.score;
}
// Setting indel, mutation, duplication/contraction
void Score::set_indel(int len, const Params &par){
	if(score < 0){cerr << "Error: set_indel for minus scores" << endl;}
	score += par.indel_val() * len;
	indel += len;
}
// Not used (for what?? I forgot..)
// void Score::set_mut(const Params &par){
// 	if(score < 0){cerr << "Error: set_mut for minus scores" << endl;}
// 	score += par.mut_val();
// 	mut++;
// }
// void Score::set_dup(int unit_id, const Params &par){
// 	if(score < 0){cerr << "Error: set_dup for minus scores" << endl;}
// 	score += par.dup_val() * len;
// 	if(dup.find(unit_id) == dup.end()){
// 		dup.emplace(unit_id, 1);
// 	}else{
// 		dup[unit_id]++;
// 	}
// }

// 従来の編集距離
Score edit_distance(const vector<int> &s, const vector<int> &t, double mut, double indel, int st, int ed){
	if(ed == -1){ed = s.size();}
	int n = ed - st;
	int m = t.size();
	
	vector<vector<Score>> dp(n+1, vector<Score>(m+1, Score(0.0)));
	for(int i = 0; i <= n; i++){dp[i][0].set_indel_2(i, indel);}
	for(int i = 0; i <= m; i++){dp[0][i].set_indel_2(i, indel);}
	
	for(int i = 1; i <= n; i++){for(int j = 1; j <= m; j++){
		double mut_ij = (s[st+i-1] == t[j-1]) ? 0.0 : mut;
		double score_ij = min({
			dp[i-1][j-1].get_score() + mut_ij,
			dp[i-1][j].get_score() + indel,
			dp[i][j-1].get_score() + indel
		});
		if(score_ij == dp[i-1][j-1].get_score() + mut_ij){
			if(mut_ij > 0.0){dp[i][j].set_mut_2(mut);}
			dp[i][j] += dp[i-1][j-1];
		}else if(score_ij == dp[i-1][j].get_score() + indel){
			dp[i][j].set_indel_2(1, indel);
			dp[i][j] += dp[i-1][j];
		}else{
			dp[i][j].set_indel_2(1, indel);
			dp[i][j] += dp[i][j-1];
		}
	}}
	
	return dp[n][m];
}

// Params class functions
// dup_pattern ... 0: 0.5*unit_length, 1: all 0.5
Params::Params(double m, double i, int dup_pattern, const vector<vector<int>> &units){
	int unit_num = units.size();
	mut = m;
	indel = i;

	// Dup score of each unit (regards dup/cont of a single unit as a indel, if not the base is included in units)
	dup_scores.resize(unit_num, Score(0.0));
	for(int i = 0; i < unit_num; i++){
		if(i < 4){
			dup_scores[i].set_indel_2(1, indel);
		}else{
			if(dup_pattern == 0) dup_scores[i].set_dup_2(i, 0.5 * units[i].size());
			else if(dup_pattern == 1) dup_scores[i].set_dup_2(i, 0.5);
			else cerr << "Error: invalid dup_pattern" << endl;
		}
	}

	// (unit, base) -> (unit, base)
	// 250424 適当に変えた。要確認。（重複の許可）
	unit_to_unit.resize(unit_num, vector<Score>(units.size()));
	for(int i = 0; i < unit_num; i++){
		for(int j = i; j < unit_num; j++){
			unit_to_unit[i][j] = min(
				edit_distance(units[i], units[j], mut, indel),
				dup_scores[i] + dup_scores[j]
			);
			unit_to_unit[j][i] = unit_to_unit[i][j];
		}
	}

	// Indel score of each unit
	indel_scores.resize(unit_num, Score(0.0));
	for(int i = 0; i < unit_num; i++){
		indel_scores[i].set_indel_2(units[i].size(), indel);
	}
}
Score Params::get_unit_to_unit(int i, int j) const {
	return unit_to_unit[i][j];
}
Score Params::get_dup(int a) const {
	return dup_scores[a];
}
Score Params::get_indel(int a) const {
	return indel_scores[a];
}
double Params::mut_val() const {return mut;}
double Params::indel_val() const {return indel;}
// double Params::dup_val() const {return dup;}

// Adds bases(A,T,G,C) and epsilon to units
void add_base_eps(
	vector<vector<int>> &units
){
	for(int i = 0; i < 4; i++){
		units.push_back({i});
	}
	units.push_back({});
}

// // debug
// void print_vec2(const vector<vector<int>> &vec2){
// 	for(auto vec : vec2){
// 		for(int val : vec){cout << val << " ";}
// 		cout << endl;
// 	}
// 	cout << endl;
// }
// void print_vec2_score(const vector<vector<Score>> &vec2){
// 	for(auto vec : vec2){
// 		for(Score val : vec){cout << val.get_score() << " ";}
// 		cout << endl;
// 	}
// 	cout << endl;
// }
// void print_vec3_score(vector<vector<vector<Score>>> vec){
// 	for(int a = 0; a < (int)vec.size(); a++){
// 		cout << "number " << a << ": " << endl;
// 		for(vector<Score> inner_vec : vec[a]){
// 			for(Score score : inner_vec){cout << setw(2) << setfill(' ') << score.get_score() << " ";}
// 			cout << endl;
// 		}
// 	}
// 	cout << endl;
// }


/* ---------------- Naive EDDC Starts from Here ---------------- */

// Below are the almost same as in "../eddc_original/eddc.h"
// Stage1 (Almost same as the original paper)
// s and t are in reverse direction, so the calculation is the same
void string_to_unit(
	const vector<int> &s,
	const vector<vector<int>> &units,
	vector<vector<Score>> &S_eps,
	vector<vector<vector<Score>>> &S,
	vector<vector<vector<Score>>> &S2,
	const Params &params
){
	int n = s.size();
	int l = units.size();
	if(n == 0){return;}
	
	// Initialization
	for(int i = 0; i < n; i++){
		S_eps[i][i+1].set_zero();
		S_eps[i][i+1].set_indel(1, params);
		for(int a = 0; a < l; a++){
			S[a][i][i+1] = params.get_unit_to_unit(a, s[i]);
		}
	}
	// for(int i = 0; i <= n; i++){
	// 	S_eps[i][i].set_zero();
	// 	for(int a = 0; a < l; a++){
	// 		S[a][i][i] = params.get_dup(a);
	// 	}
	// }
	
	// DP
	for(int j = 2; j <= n; j++){
		for(int i = j-2; i >= 0; i--){
			// Update S2
			for(int a = 0; a < l; a++){
				for(int h = i+1; h < j; h++){
					Score tmp = min({
						S[a][i][h] + S_eps[h][j],
						S_eps[i][h] + S[a][h][j],
						params.get_dup(a) + S[a][i][h] + S[a][h][j]
					});

					// Added the part below, for the original algorithm does not work well
					if(j-i <= (int)units[a].size()){
						tmp = min(tmp, edit_distance(s, units[a], params.mut_val(), params.indel_val(), i, j));
					}

					if(S2[a][i][j].get_score() < 0){
						S2[a][i][j] = tmp;
					}else{
						S2[a][i][j] = min(S2[a][i][j], tmp);
					}
				}
			}
			// Update S
			for(int a = 0; a < l; a++){
				for(int b = 0; b < l; b++){
					if(S[a][i][j].get_score() < 0){
						S[a][i][j] = params.get_unit_to_unit(a,b) + S2[b][i][j];
					}else{
						S[a][i][j] = min(S[a][i][j], params.get_unit_to_unit(a,b) + S2[b][i][j]);
					}
				}
			}
			// Update S_eps
			for(int a = 0; a < l; a++){
				Score cand = min(params.get_dup(a) + S[a][i][j], params.get_indel(a) + S[a][i][j]);
				if(S_eps[i][j].get_score() < 0){
					S_eps[i][j] = cand;
				}else{
					S_eps[i][j] = min(S_eps[i][j], cand);
				}
			}
		}
	}
}

Score calc_eddc(
	const vector<int> &s,
	const vector<int> &t,
	const vector<vector<int>> &units,
	const Params &params
){
	int n = s.size();
	int m = t.size();
	int l = units.size();
	if(n == 0 && m == 0){return Score(0.0);}
	
	vector<vector<Score>> ED(n+1, vector<Score>(m+1, Score(-1.0)));
	vector<vector<vector<Score>>> EDT(l, vector<vector<Score>>(n+1, vector<Score>(m+1, Score(-1.0))));
	
	vector<vector<Score>> S_eps(n+1, vector<Score>(n+1, Score(-1.0)));
	vector<vector<vector<Score>>> S(l, vector<vector<Score>>(n+1, vector<Score>(n+1, Score(-1.0))));
	vector<vector<vector<Score>>> S2(l, vector<vector<Score>>(n+1, vector<Score>(n+1, Score(-1.0))));
	vector<vector<Score>> T_eps(m+1, vector<Score>(m+1, Score(-1.0)));
	vector<vector<vector<Score>>> T(l, vector<vector<Score>>(m+1, vector<Score>(m+1, Score(-1.0))));
	vector<vector<vector<Score>>> T2(l, vector<vector<Score>>(m+1, vector<Score>(m+1, Score(-1.0))));
	
	// Stage1
	string_to_unit(s, units, S_eps, S, S2, params);
	// if(n == 3 && s[1] == 3 && t[0] == 1){print_vec3_score(S);}
	// print_vec2(S_eps);
	// print_vec3(S);
	// print_vec3(S2);
	string_to_unit(t, units, T_eps, T, T2, params);
	// cout << "T_eps "; print_vec2(T_eps);
	// cout << "T "; print_vec3(T);
	// cout << "T2 "; print_vec3(T2);

	if(n == 0){return T_eps[0][m];}
	if(m == 0){return S_eps[0][n];}

	// Stage2
	// Initialization
	ED[0][0].set_zero();
	for(int i = 1; i < m+1; i++){
		ED[0][i] = T_eps[0][i];
		ED[1][i] = T[s[0]][0][i];
	}
	for(int i = 1; i < n+1; i++){
		ED[i][0] = S_eps[0][i];
		ED[i][1] = S[t[0]][0][i];
	}

	// Initialization of EDT (i=1)
	for(int a = 0; a < l; a++){
		for(int j = 2; j < m+1; j++){
			for(int h = 1; h < j; h++){
				if(EDT[a][1][j].get_score() < 0){EDT[a][1][j] = ED[1][h] + T[a][h][j];}
				else{EDT[a][1][j] = min(EDT[a][1][j], ED[1][h] + T[a][h][j]);}
			}
		}
	}

	// DP
	for(int j = 2; j <= m; j++){
		for(int i = 2; i <= n; i++){
			// Update EDT
			for(int a = 0; a < l; a++){
				for(int h = 1; h < j; h++){
					if(EDT[a][i][j].get_score() < 0){
						EDT[a][i][j] = ED[i][h] + T[a][h][j];
					}else{
						EDT[a][i][j] = min(EDT[a][i][j], ED[i][h] + T[a][h][j]);
					}
					
					// cout << "EDT[a][i][j]: " << EDT[a][i][j].score << endl;
					// cout << "ED[i][h]: " << ED[i][h].score << endl;
					// cout << "T[a][h][j]: " << T[a][h][j].score << endl;
				}
			}
			// Update ED
			for(int h = 1; h < i; h++){
				for(int a = 0; a < l; a++){
					Score tmp = min(S[a][0][i] + T[a][0][j], EDT[a][h][j] + S[a][h][i]);
					if(ED[i][j].get_score() < 0){
						ED[i][j] = tmp;
					}else{
						ED[i][j] = min(ED[i][j], tmp);
					}
					
					// cout << "(i,j,h,a)=(" << i << "," << j << "," << h << "," << a << ")" << endl;
					// cout << "S[a][0][i]: " << S[a][0][i].score << endl;
					// cout << "T[a][0][j]: " << T[a][0][j].score << endl;
					// cout << "EDT[a][h][j]: " << EDT[a][h][j].score << endl;
					// cout << "S[a][h][i]: " << S[a][h][i].score << endl;
					// cout << endl;
				}
			}
		}
	}
	
	// cout << "ED "; print_vec2(ED);
	// cout << "EDT "; print_vec3(EDT);
	
	return ED[n][m];
}

/* ---------------- Naive EDDC Ends Here ---------------- */


// Creates eddc_units（仮）（あとで共通化する）
void shift_units(const vector<vector<int>> &units, vector<vector<int>> &eddc_units){
	int ub_num = units.size();
	int u_num = ub_num - 5;

	eddc_units.clear();
	eddc_units.reserve(ub_num - 1);

	// bases
	for(int i = 0; i < 4; i++){
		eddc_units.push_back({i});
	}
	// units
	for(int i = 0; i < u_num; i++){
		eddc_units.push_back(units[i]);
	}
}

void calc_f(
	const vector<vector<int>> &encodings,
	const vector<vector<int>> &units,
	const Params &params,
	vector<vector<vector<vector<Score>>>> &f_scores,
	vector<int> acc_par1s,
	vector<long long> &measure_time
){
	int ub_num = units.size();
	int reads_num = encodings.size();
	int eps = ub_num - 1;

	// calc_eddc()の仕様上, eddc_unitsは[{bases}, {units}]の順番にする必要がある。これはあとで修正する。
	vector<vector<int>> eddc_units;
	shift_units(units, eddc_units);

	auto c1c2_st = chrono::system_clock::now();

	// c1[x,y]		: calc_eddc(x,y)
	vector<vector<Score>> c1(ub_num, vector<Score>(ub_num, Score(0.0)));
	for(int x = 0; x < ub_num; x++){
		for(int y = x+1; y < ub_num; y++){
			c1[x][y] = calc_eddc(units[x], units[y], eddc_units, params);
			c1[y][x] = c1[x][y];
		}
	}

	// c2[x,y,z]	: calc_eddc(xy,z)
	vector<vector<vector<Score>>> c2(
		ub_num, vector<vector<Score>>(
			ub_num, vector<Score>(ub_num, Score(DBL_MAX)))
	);
	for(int x = 0; x < ub_num; x++){
		for(int y = 0; y < ub_num; y++){
			for(int z = 0; z < ub_num; z++){
				vector<int> xy;
				xy.reserve(units[x].size() + units[y].size());
				xy.insert(xy.end(), units[x].begin(), units[x].end());
				xy.insert(xy.end(), units[y].begin(), units[y].end());
				c2[x][y][z] = calc_eddc(xy, units[z], eddc_units, params);
			}
		}
	}

	auto c1c2_ed = chrono::system_clock::now();
	auto c1c2 = chrono::duration_cast<chrono::milliseconds>(c1c2_ed - c1c2_st).count();
	if(MEASURE_TIME){cout << "calc c1, c2: " << c1c2 << " msec" << endl;}

	// cout << "c2: " << endl;
	// print_vec3_score(c2);

	// valid_rules	: xy -> z st. \forall {p,q} {eddc(xy, z) <= eddc(x,p) + eddc(y,q) + eddc(pq,z)}
	// 各valid_rule xy -> z に対して、<x,z>のsetと<x,y,z>のset
	set<pair<int,int>> valid_rules_xz;
	set<tuple<int,int,int>> valid_rules_xyz;
	for(int x = 0; x < ub_num; x++){
		for(int y = 0; y < ub_num; y++){
			for(int z = 0; z < ub_num; z++){
				bool is_redundant = false;
				for(int p = 0; p < ub_num; p++){
					for(int q = 0; q < ub_num; q++){
						if(p == x && q == y){continue;}
						if(c2[x][y][z] >= c1[x][p] + c1[y][q] + c2[p][q][z]){
							is_redundant = true;
							break;
						}
					}
					if(is_redundant){break;}
				}
				if(!is_redundant){
					valid_rules_xz.insert(make_pair(x,z));
					valid_rules_xyz.insert(make_tuple(x,y,z));
				}
			}
		}
	}

	auto vrule_ed = chrono::system_clock::now();
	auto vrule = chrono::duration_cast<chrono::milliseconds>(vrule_ed - c1c2_ed).count();
	if(MEASURE_TIME){cout << "valid rules: " << vrule << " msec" << endl;}
	measure_time[0] = chrono::duration_cast<chrono::milliseconds>(vrule_ed - c1c2_st).count();

	// cout << "calc valid" << endl;
	// for(auto &[x,y,z] : valid_rules_xyz){cout << x << y << z << endl;}

	// f[w,b,e,z]
	for(int w = 0; w < reads_num; w++){
		int n = encodings[w].size();

		// dp[b,e,z]: f_scores[w][b,e,z] = Cost(w[b,e), z) = min_{b<k<=e}{min_{x\in U}{dp[b,k,x] + dp_sub[k,e,x,z]}}
		// eをlen (= e-b) に変更
		vector<vector<vector<Score>>> dp(
			n+1, vector<vector<Score>>(
				acc_par1s[w]+1, vector<Score>(ub_num, Score(DBL_MAX)))
		);
		// dp_sub[k,e,x,z] = min_{y\in U}{f[w,k,e,y] + eddc(xy,z)}
		// eをlen (= e-k) に変更
		vector<vector<vector<vector<Score>>>> dp_sub(
			n+1, vector<vector<vector<Score>>>(
				acc_par1s[w]+1, vector<vector<Score>>(
					ub_num, vector<Score>(ub_num, Score(DBL_MAX))))
		);

		// Initialization of dp
		for(int b = 0; b < n+1; b++){
			// b == eのとき, dp[b,e,z] = c1[eps,z]
			for(int z = 0; z < ub_num; z++){dp[b][0][z] = c1[eps][z];}
			// b+1 == eのとき, dp[b,e,z] = c1[b[b,e),z] = c1[b[b],z]
			if(b != n){
				for(int z = 0; z < ub_num; z++){
					dp[b][1][z] = c1[encodings[w][b]][z];
				}
			}
		}

		// Initialization of dp_sub
		for(int b = 0; b < n+1; b++){
			for(int len = 0; len < 2 && b+len < n+1; len++){
				for(int z = 0; z < ub_num; z++){
					for(int x = 0; x < ub_num; x++){
						for(int y = 0; y < ub_num; y++){
							dp_sub[b][len][x][z] = min(dp_sub[b][len][x][z], dp[b][len][y] + c2[x][y][z]);
						}
					}
				}
			}
		}

		// DP (CYK-like algorithm)
		for(int len = 2; len < acc_par1s[w]+1; len++){
			for(int b = 0; b < n+1-len; b++){
				// dp[b,e,z] = min_{b<k<=e}{min_{x\in U}{dp[b,k,x] + dp_sub[k,e,x,z]}}
				// kをbからのindexに変更
				for(int k = 1; k <= len; k++){
					for(const auto &[x,z] : valid_rules_xz){
						dp[b][len][z] = min(dp[b][len][z], dp[b][k][x] + dp_sub[b+k][len-k][x][z]);
					}
				}
				// dp_sub[k,e,x,z] = min_{y\in U}{dp[k,e,y] + eddc(xy,z)}
				for(const auto &[x,y,z] : valid_rules_xyz){
					dp_sub[b][len][x][z] = min(dp_sub[b][len][x][z], dp[b][len][y] + c2[x][y][z]);
				}
			}
		}

		f_scores[w] = dp;
	}

	auto f_dp_ed = chrono::system_clock::now();
	auto f_dp = chrono::duration_cast<chrono::milliseconds>(f_dp_ed - vrule_ed).count();
	if(MEASURE_TIME){cout << "f main dp: " << f_dp << " msec" << endl;}
	measure_time[1] = f_dp;
}

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
){
	int n = u.size();
	int m = v.size();
	int ub_num = units.size();

	// dp[i,j] = Cost(u[0,i), v[0,j))
	vector<vector<Score>> dp(n+1, vector<Score>(m+1, Score(DBL_MAX)));
	// dp_sub[p,j,x] = min_{0<=q<j}{dp[p,q] + f[v,q,j,x]}
	vector<vector<vector<Score>>> dp_sub(n+1, vector<vector<Score>>(m+1, vector<Score>(ub_num, Score(DBL_MAX))));

	// Initialization
	dp[0][0] = Score(0.0);

	// DP
	for(int i = 0; i < n+1; i++){
		for(int j = 0; j < m+1; j++){
			if(abs(i - j) > acc_par2){continue;}
			// dp_sub[i,j,x]
			for(int x = 0; x < ub_num; x++){
				Score val(DBL_MAX);
				// dp_sub[p,j,x] = min_{0<=q<j}{dp[p,q] + f[v,q,j,x]}	(p=i)
				for(int q = max(0, j-acc_par1s[v_idx]); q <= j; q++){	// とにかくf[-,q,j,-]で[q,j)の長さをacc_par1で制限すればよい
					if(abs(i - q) > acc_par2){continue;}
					// i==p && j==q のときはDBL_MAXのままなのでスキップできる
					val = min(val, dp[i][q] + f_scores[v_idx][q][j-q][x]);
				}
				dp_sub[i][j][x] = val;
			}

			// dp[i,j] = min_{x\in U}{min_{0<=p<i}{f[u,p,i,x] + dp_sub[p,j,x]}}
			Score val(DBL_MAX);
			for(int x = 0; x < ub_num; x++){
				for(int p = max(0, i-acc_par1s[u_idx]); p <= i; p++){
					if(abs(p - j) > acc_par2){continue;}
					val = min(val, f_scores[u_idx][p][i-p][x] + dp_sub[p][j][x]);
				}
			}
			dp[i][j] = val;
		}
	}

	return dp[n][m];
}

// 全てのTRの組に対してEDDCを計算する
void heddc_all(
	const vector<vector<int>> &encodings,
	vector<vector<int>> &units,
	const Params &params,
	vector<vector<Score>> &scores,
	int par,
	vector<long long> &measure_time 	// c1 c2 v_rules, f, maind dp の実行時間
){
	int tr_num = encodings.size();
	measure_time.resize(3);

	// f_scores計算時の最大長さ
	// 最小値を設定
	int acc_par1 = par;
	// acc1を [encoding長 / 最短encoding長] を基準にする
	int min_enc_len = INT_MAX;
	for(const auto &enc : encodings){
		if((int)enc.size() < min_enc_len){min_enc_len = enc.size();}
	}
	min_enc_len = max(1, min_enc_len-1);    // バッファ (& 0にならないように)
	vector<int> acc_par1s(tr_num);
	for(int i = 0; i < tr_num; i++){
		acc_par1s[i] = max(acc_par1, int(encodings[i].size() / min_enc_len));
		acc_par1s[i] = max(acc_par1s[i], (int)encodings[i].size() - min_enc_len);
	}
	// // main dpでの最大長さ (が[read長の差 / 最小ユニット長]の何倍か)
	// // enc長の差の何倍かに変更
	// double acc_par2_ratio = (double)INT_MAX / 10.0;
	// int min_unit_len = INT_MAX;
	// for(const auto &unit : units){
	// 	if(unit.size() < min_unit_len){min_unit_len = unit.size();}
	// }

	// unitsにbaseとepsを追加する
	add_base_eps(units);
	// int ub_num = units.size();

	// f[u,b,e,x] = Pr(u[b,e) -> x)
	// 各u (read) に対してfを計算する
	vector<vector<vector<vector<Score>>>> f_scores(tr_num);
	calc_f(encodings, units, params, f_scores, acc_par1s, measure_time);

	// 各組み合わせについてheddcを実行する
	scores.assign(tr_num, vector<Score>(tr_num, Score(0.0)));
	auto dp_st = chrono::system_clock::now();
	for(int i = 0; i < tr_num; i++){
		for(int j = i+1; j < tr_num; j++){
			// int acc_par2 = int(
			// 	abs((int)encodings[i].size() - (int)encodings[j].size()) * acc_par2_ratio
			// );
			// acc_par2 = max(acc_par2, int((double)max(encodings[i].size(), encodings[j].size())/5));
			// double val = calc_heddc(encodings[i], i, encodings[j], j, units, params, f_scores, acc_par1s, acc_par2);
			Score val = calc_heddc(encodings[i], i, encodings[j], j, units, params, f_scores, acc_par1s, INT_MAX);
			scores[i][j] = val;
			scores[j][i] = scores[i][j];
		}
	}
	auto dp_ed = chrono::system_clock::now();
	auto heddc_all_dp = chrono::duration_cast<chrono::milliseconds>(dp_ed - dp_st).count();
	if(MEASURE_TIME){cout << "heddc_all main dp: " << heddc_all_dp << " msec" << endl;}
	measure_time[2] = heddc_all_dp;
}
