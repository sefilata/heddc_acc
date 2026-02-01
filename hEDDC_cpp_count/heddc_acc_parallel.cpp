#include "heddc_acc_parallel.h"

// Constructor of Score
Score::Score() : score(0), mut(0), indel(0), dup() {}
Score::Score(double s) : score(s), mut(0), indel(0), dup() {}
// Setting zero for Score
void Score::set_zero(){score = 0; mut = 0; indel = 0; dup.clear();}
// Getting variables in Score
double Score::get_score() const {return score;}
int Score::get_mut() const {return mut;}
int Score::get_indel() const {return indel;}
void Score::print_dup(ostream& os, const vector<vector<int>> &units) const {
	for(const auto &dup_ : dup){
		if(dup_.second > 0){
			if(dup_.first < 4){
				// cerr << "Warning: print_dup() called with unit_id < 4, which is not expected." << endl;
				os << "(";
				os << base_map_reverse[dup_.first];
				os << ", " << dup_.second << ")";
				continue;
			}
			os << "(";
			string unit_str;
			digits_to_string(units[dup_.first], unit_str);
			os << unit_str;
			os << ", " << dup_.second << ")";
		}
	}
}
// Used only in edit_distance() and Constructor of Params
void Score::set_indel_2(int len, double par){score += par * len; indel += len;}
void Score::set_mut_2(double par){score += par; mut++;}
void Score::set_dup_2(int unit_id, double par){score += par;; dup.emplace(unit_id, 1);}
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
bool Score::operator==(const Score& other) const {
	return fabs(score - other.score) < 1e-12;
}
// Setting indel, mutation, duplication/contraction
void Score::set_indel(int len, const Params &par){
	if(score < 0){cerr << "Error: set_indel for minus scores" << endl;}
	score += par.indel_val() * len;
	indel += len;
}

// Conventional edit distance: ed(s[st,ed), t)
Score edit_distance(const vector<int> &s, const vector<int> &t, double mut, double indel, int st, int ed){
	if(ed == -1){ed = s.size();}
	int n = ed - st;
	int m = t.size();
	
	vector<vector<Score>> dp(n+1, vector<Score>(m+1, Score(0.0)));
	for(int i = 0; i <= n; i++){dp[i][0].set_indel_2(i, indel);}
	for(int i = 0; i <= m; i++){dp[0][i].set_indel_2(i, indel);}
	
	for(int i = 1; i <= n; i++){for(int j = 1; j <= m; j++){
		double mut_ij =(s[st+i-1] == t[j-1]) ? 0.0 : mut;
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
// dup_pattern ... 0: PARAM_0_VALUE*unit_length, 1: all PARAM_1_VALUE
// units: [{bases}. {units}, eps]
Params::Params(double m, double i, int dup_pattern, const vector<vector<int>> &units){
	int unit_num = units.size();
	mut = m;
	indel = i;

	// Dup score of each unit (regarding dup/cont of a single unit as indel, if the base is not included as units)
	dup_scores.resize(unit_num, Score(0.0));
	for(int i = 0; i < unit_num; i++){
		if(i < 4){
			dup_scores[i].set_indel_2(1, indel);
		}else{
			if(dup_pattern == 0) dup_scores[i].set_dup_2(i, PARAM_0_VALUE * units[i].size());
			else if(dup_pattern == 1) dup_scores[i].set_dup_2(i, PARAM_1_VALUE);
			else cerr << "Error: invalid dup_pattern" << endl;
		}
	}

	// ed(x, y): これいらないかも？最後に要確認
	unit_to_unit.resize(unit_num, vector<Score>(unit_num));
	for(int i = 0; i < unit_num; i++){
		for(int j = i; j < unit_num; j++){
			unit_to_unit[i][j] = edit_distance(units[i], units[j], mut, indel);
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

// DEBUG
void print_vec2_score(const vector<vector<Score>> &vec2, const vector<vector<int>> &units){
	for(auto vec : vec2){
		for(Score val : vec){cout << val.get_score() << " ";}
		cout << endl;
	}
	cout << "variants" << endl;
	for(auto vec : vec2){
		for(Score val : vec){
			cout << "{mut:" << val.get_mut();
			cout << ", indel:" << val.get_indel();
			cout << ", dup:";
			val.print_dup(cout, units);
			cout << "} ";
		}
		cout << endl;
	}
	cout << endl;
}

// Measure memory (Linux only)
long long get_hwm_kb(){
	ifstream f("/proc/self/status");
	string line;
	while(getline(f, line)){
		if(line.rfind("VmHWM:", 0) == 0){
			long long v; string unit;
			istringstream iss(line.substr(6));
			iss >> v >> unit;
			return v; 	// kB
		}
	}
	return -1; 	// error
}


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
){
	int n = s.size();
	int l = units.size();
	// int eps = l - 1;
	if(n == 0){return;}
	
	// Initialization
	for(int i = 0; i < n; i++){
		S_eps[i][i+1].set_zero();
		S_eps[i][i+1].set_indel(1, params);
		for(int a = 0; a < l; a++){
			S[a][i][i+1] = params.get_unit_to_unit(a, s[i]);
		}
	}
	
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
					// I'm not sure how long it should be computed..
					// if(j-i <= (int)units[a].size()){
					// 	tmp = min(tmp, edit_distance(s, units[a], params.mut_val(), params.indel_val(), i, j));
					// }
					tmp = min(tmp, edit_distance(s, units[a], params.mut_val(), params.indel_val(), i, j));

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
				Score cand = params.get_indel(a) + S[a][i][j];
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
	string_to_unit(t, units, T_eps, T, T2, params);

	if(n == 0){return T_eps[0][m];}
	if(m == 0){return S_eps[0][n];}

	// Stage2
	// Initialization of ED
	ED[0][0].set_zero();
	for(int i = 1; i < m+1; i++){
		ED[0][i] = T_eps[0][i];
		ED[1][i] = T[s[0]][0][i];
	}
	for(int i = 1; i < n+1; i++){
		ED[i][0] = S_eps[0][i];
		ED[i][1] = S[t[0]][0][i];
	}

	// Initialization of EDT(i=1)
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
				}
			}
		}
	}
	
	return ED[n][m];
}

/* ---------------- Naive EDDC Ends Here ---------------- */


// heddc(s[i,j), x) for any combination of s\in{Encoded TR},x\in{Units},i,j
void calc_f(
	const vector<vector<int>> &encodings,
	const vector<vector<int>> &units,
	const Params &params,
	vector<vector<vector<vector<Score>>>> &f_scores,
	const vector<int> &acc_par1s,
	vector<long long> &measure_time,
	vector<long long> &measure_mem
){
	int ub_num = units.size();
	int reads_num = encodings.size();
	int eps = ub_num - 1;

	auto c1c2_st = chrono::system_clock::now();

	// Rank the lengths of units for load leveling
	vector<int> lengths_rank(ub_num);
	iota(lengths_rank.begin(), lengths_rank.end(), 0);
	sort(lengths_rank.begin(), lengths_rank.end(), [&](int a, int b){
		return units[a].size() > units[b].size();
	});

	// c1[x,y] 	: calc_eddc(x,y)
	if(DEBUG_OUT){
		auto now = chrono::system_clock::now();
		time_t t = chrono::system_clock::to_time_t(now);
		cout << "Started calc-c1 at " << ctime(&t);
		cout.flush();
	}
	vector<vector<Score>> c1(ub_num, vector<Score>(ub_num, Score(0.0)));
	#pragma omp parallel for schedule(dynamic)
	for(int xy = 0; xy < ub_num*ub_num; xy++){
		int x = lengths_rank[xy / ub_num];
		int y = lengths_rank[xy % ub_num];
		if(x < 4 || y < 4 || x >= y){continue;} 	// Compute only for units and eps (excluding bases).
		c1[x][y] = calc_eddc(units[x], units[y], units, params);
		c1[y][x] = c1[x][y];
	}

	// c2[x,y,z] 	: calc_eddc(xy,z)
	if(DEBUG_OUT){
		auto now = chrono::system_clock::now();
		time_t t = chrono::system_clock::to_time_t(now);
		cout << "Started calc-c2 at " << ctime(&t);
		cout.flush();
	}
	vector<vector<vector<Score>>> c2(
		ub_num, vector<vector<Score>>(
			ub_num, vector<Score>(ub_num, Score(DBL_MAX)))
	);
	#pragma omp parallel for schedule(dynamic)
	for(int xyz = 0; xyz < ub_num*ub_num*ub_num; xyz++){
		int x = lengths_rank[xyz /(ub_num * ub_num)];
		int y = lengths_rank[(xyz / ub_num) % ub_num];
		int z = lengths_rank[xyz % ub_num];
		if(x < 4 || y < 4 || z < 4){continue;}
		// 2025-11-27 	In the case x==y==z, omit the computation.（Warn: It could be inappropriate for certain parameter settings.）
		if(x == y && y == z){
			c2[x][y][z] = params.get_dup(x);
			continue;
		}
		vector<int> xy;
		xy.reserve(units[x].size() + units[y].size());
		xy.insert(xy.end(), units[x].begin(), units[x].end());
		xy.insert(xy.end(), units[y].begin(), units[y].end());
		c2[x][y][z] = calc_eddc(xy, units[z], units, params);
	}

	auto c1c2_ed = chrono::system_clock::now();
	auto c1c2 = chrono::duration_cast<chrono::milliseconds>(c1c2_ed - c1c2_st).count();
	if(MEASURE_TIME){cout << "calc c1, c2: " << c1c2 << " msec" << endl;}

	// valid_rules	: xy -> z st. \forall {p,q} {eddc(xy, z) <= eddc(x,p) + eddc(y,q) + eddc(pq,z)}
	// For each valid_rule xy -> z, store the set of <x,z> and <x,y,z>
	if(DEBUG_OUT){
		auto now = chrono::system_clock::now();
		time_t t = chrono::system_clock::to_time_t(now);
		cout << "Started valid-rules at " << ctime(&t);
		cout.flush();
	}
	set<pair<int,int>> valid_rules_xz;
	set<tuple<int,int,int>> valid_rules_xyz;
	for(int x = 0; x < ub_num; x++){
		if(x < 4){continue;}
		for(int y = 0; y < ub_num; y++){
			if(y < 4){continue;}
			for(int z = 0; z < ub_num; z++){
				if(z < 4){continue;}
				bool is_redundant = false;
				for(int p = 0; p < ub_num; p++){
					if(p < 4){continue;}
					for(int q = 0; q < ub_num; q++){
						if(q < 4){continue;}
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
	if(MEASURE_MEMORY_LINUX) measure_mem[0] = get_hwm_kb();

	// f[w,b,e,z] = heddc(w[b,e),z)
	if(DEBUG_OUT){
		auto now = chrono::system_clock::now();
		time_t t = chrono::system_clock::to_time_t(now);
		cout << "Started f-scores at " << ctime(&t);
		cout.flush();
	}
	#pragma omp parallel for schedule(dynamic)
	for(int w = 0; w < reads_num; w++){
		int n = encodings[w].size();
		int tab_size = min(n, acc_par1s[w])+1;

		// dp[b,e,z]: f_scores[w][b,e,z] = Cost(w[b,e), z) = min_{b<k<=e}{min_{x\in U}{dp[b,k,x] + dp_sub[k,e,x,z]}}
		// eをlen(= e-b) に変更
		vector<vector<vector<Score>>> dp(
			n+1, vector<vector<Score>>(
				tab_size, vector<Score>(ub_num, Score(DBL_MAX)))
		);
		// dp_sub[k,e,x,z] = min_{y\in U}{f[w,k,e,y] + eddc(xy,z)}
		// eをlen(= e-k) に変更
		vector<vector<vector<vector<Score>>>> dp_sub(
			n+1, vector<vector<vector<Score>>>(
				tab_size, vector<vector<Score>>(
					ub_num, vector<Score>(ub_num, Score(DBL_MAX))))
		);

		// Initialization of dp
		for(int b = 0; b < n+1; b++){
			// Case b == e, dp[b,e,z] = c1[eps,z]
			for(int z = 0; z < ub_num; z++){
				if(z < 4){continue;}
				dp[b][0][z] = c1[eps][z];
			}
			// Case b+1 == e, dp[b,e,z] = c1[b[b,e),z] = c1[b[b],z]
			if(b != n){
				for(int z = 0; z < ub_num; z++){
					if(z < 4){continue;}
					dp[b][1][z] = c1[encodings[w][b]][z];
				}
			}
		}

		// Initialization of dp_sub
		for(int b = 0; b < n+1; b++){
			for(int len = 0; len < 2 && b+len < n+1; len++){
				for(int z = 0; z < ub_num; z++){
					if(z < 4){continue;}
					for(int x = 0; x < ub_num; x++){
						if(x < 4){continue;}
						for(int y = 0; y < ub_num; y++){
							if(y < 4){continue;}
							dp_sub[b][len][x][z] = min(dp_sub[b][len][x][z], dp[b][len][y] + c2[x][y][z]);
						}
					}
				}
			}
		}

		// DP(CYK-like algorithm)
		for(int len = 2; len < min(n, acc_par1s[w])+1; len++){
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
	if(MEASURE_MEMORY_LINUX) measure_mem[1] = get_hwm_kb();
	measure_time[1] = f_dp;
}

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
){
	int n = u.size();
	int m = v.size();
	int ub_num = units.size();
	// int eps = ub_num - 1;

	// dp[i,j] = heddc(u[0,i), v[0,j))
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
				if(x < 4){continue;}
				Score val(DBL_MAX);
				// dp_sub[p,j,x] = min_{0<=q<j}{dp[p,q] + f[v,q,j,x]}	(p=i)
				for(int q = max(0, j-acc_par1s[v_idx]); q <= j; q++){ 	// Limit the length of [q,j) in f[-,q,j,-] by acc_par1
					if(abs(i - q) > acc_par2){continue;}
					// In case i==p && j==q, it is skipped since val remains to be DBL_MAX
					val = min(val, dp[i][q] + f_scores[v_idx][q][j-q][x]);
				}
				dp_sub[i][j][x] = val;
			}

			// dp[i,j] = min_{x\in U}{min_{0<=p<i}{f[u,p,i,x] + dp_sub[p,j,x]}}
			Score val(DBL_MAX);
			for(int x = 0; x < ub_num; x++){
				if(x < 4){continue;}
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

// Compute hEDDC for all TR pairs
void heddc_all(
	const vector<vector<int>> &encodings,
	vector<vector<int>> &units,
	const Params &params,
	vector<vector<Score>> &scores,
	int par,
	vector<long long> &measure_time, 	// Computation times of c1_c2 v_rules, f, main_dp
	vector<long long> &measure_mem, 	// Max memory consumption when c1_c2, f, or main_dp finished
	int parallel_num_1, 	// # of threads other than main dp
	int parallel_num_2, 	// # of threads in main dp
	vector<int> &tr_nums 	// # of sequences before and after compression (consolidation of seqs with identical unit decompositions)
){
	int tr_num = encodings.size();
	measure_time.resize(3);
	omp_set_num_threads(parallel_num_1);

	// skip calculating for same encodings
	// calculating index when compressed
	vector<int> index_map(tr_num, -1); 	// index when compressed
	int tr_num_compressed = 0;
	for(int i = 0; i < tr_num; i++){
		bool found = false;
		for(int j = 0; j < i; j++){
			if(encodings[i] != encodings[j]) continue;
			index_map[i] = index_map[j];
			found = true;
			break;
		}
		if(found) continue;
		index_map[i] = tr_num_compressed;
		tr_num_compressed++;
	}
	if(DEBUG_OUT) cout << "[INFO] tr_num: " << tr_num << endl;
	if(DEBUG_OUT) cout << "[INFO] tr_num_compressed: " << tr_num_compressed << endl;
	// create compressed encodings vector
	vector<vector<int>> encodings_compressed(tr_num_compressed);
	vector<bool> filled(tr_num_compressed, false);
	for(int i = 0; i < tr_num; i++){
		if(filled[index_map[i]]) continue;
		encodings_compressed[index_map[i]] = encodings[i];
		filled[index_map[i]] = true;
	}

	// Length limits for f_scores computations
	// Setting the minimum value of length limits
	vector<int> acc_par1s(tr_num_compressed, par);
	// par == INT_MAX if ACC_PATTERN is 0
	if(par != INT_MAX){
		// (i) TRs containing corresponding repeats with largely different lengths.
		// acc1 is set to [{length of target encoding} / {length of shortest encoding}]
		// (ii) A TR contains units that are absent in the other TR.
		// acc1 is set to [{length of target encoding} - {length of shortest encoding}]
		int min_enc_len = INT_MAX;
		for(const auto &enc : encodings_compressed){
			if((int)enc.size() < min_enc_len){min_enc_len = enc.size();}
		}
		min_enc_len = max(1, min_enc_len-1);	// Buffer & avoiding 0-division
		for(int i = 0; i < tr_num_compressed; i++){
			acc_par1s[i] = max(acc_par1s[i], int(encodings_compressed[i].size() / min_enc_len));
			acc_par1s[i] = max(acc_par1s[i],(int)encodings_compressed[i].size() - min_enc_len);
		}
	}

	// f[u,b,e,x] = heddc(u[b,e), x)
	// Compute f for each u \in{encoded TRs}
	vector<vector<vector<vector<Score>>>> f_scores(tr_num);
	calc_f(encodings_compressed, units, params, f_scores, acc_par1s, measure_time, measure_mem);

	// Rank the lengths of units for leveling the loads
	vector<int> lengths_rank(tr_num_compressed);
	iota(lengths_rank.begin(), lengths_rank.end(), 0);
	sort(lengths_rank.begin(), lengths_rank.end(), [&](int a, int b){
		return encodings_compressed[a].size() > encodings_compressed[b].size();
	});

	// Compute main_dp of heddc
	omp_set_num_threads(parallel_num_2);
	if(DEBUG_OUT){
		auto now = chrono::system_clock::now();
		time_t t = chrono::system_clock::to_time_t(now);
		cout << "Started main-dp at " << ctime(&t);
		cout.flush();
	}
	vector<vector<Score>> scores_compressed(tr_num_compressed);
	scores_compressed.assign(tr_num_compressed, vector<Score>(tr_num_compressed, Score(0.0)));
	auto dp_st = chrono::system_clock::now();
	#pragma omp parallel for schedule(dynamic)
	for(int ij = 0; ij < tr_num_compressed * tr_num_compressed; ij++){
		int i = lengths_rank[ij / tr_num_compressed];
		int j = lengths_rank[ij % tr_num_compressed];
		if(i >= j){continue;}

		Score val = calc_heddc(encodings_compressed[i], i, encodings_compressed[j], j, units, params, f_scores, acc_par1s, INT_MAX);

		scores_compressed[i][j] = val;
		scores_compressed[j][i] = scores_compressed[i][j];
	}

	// Scores
	scores.assign(tr_num, vector<Score>(tr_num, Score(0.0)));
	#pragma omp parallel for schedule(dynamic)
	for(int ij = 0; ij < tr_num * tr_num; ij++){
		int i = ij / tr_num;
		int j = ij % tr_num;
		scores[i][j] = scores_compressed[index_map[i]][index_map[j]];
	}

	auto dp_ed = chrono::system_clock::now();
	auto heddc_all_dp = chrono::duration_cast<chrono::milliseconds>(dp_ed - dp_st).count();
	if(MEASURE_TIME){cout << "heddc_all main dp: " << heddc_all_dp << " msec" << endl;}
	if(MEASURE_MEMORY_LINUX) measure_mem[2] = get_hwm_kb();
	measure_time[2] = heddc_all_dp;
	tr_nums = {tr_num, tr_num_compressed};
}
 