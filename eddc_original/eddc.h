#pragma once
#include <iostream>
#include <vector>
#include <chrono>
#include <map>
#include <algorithm>
#include <iomanip>
using namespace std;

//#define DEBUG

class Params;

// スコア構造体
class Score{
private:
	double score;
	int mut, indel;
	map<int, int> dup;

public:
	Score() : score(0), mut(0), indel(0), dup() {}
	Score(double s) : score(s), mut(0), indel(0), dup() {}

	// 各変数の追加
	void set_indel(int len, const Params &par);
	void set_mut(const Params &par);
	void set_dup(int len, const Params &par);
	void set_zero(){score = 0;}

	// 変数の取得
	double get_score(){return score;}
	int get_mut(){return mut;}
	int get_indel(){return indel;}
	void print_dup(){
		for(auto it = dup.begin(); it != dup.end(); it++){
			if(it != dup.begin()){cout << ", ";}
			cout << "(" << it->first << ": " << it->second << ")";
		}
	}

	// edit_distance()とParamsのコンストラクタでのみ使用
	void set_indel_2(int len, double par){score += par * len; indel += len;}
	void set_mut_2(double par){score += par; mut++;}
	void set_dup_2(int len, double par){score += par * len; dup.emplace(len, 1);}

	// Score同士の加算
	Score operator+(const Score &other){
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

	// Score同士の加算代入
	Score& operator+=(const Score &other) {
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

	// Score同士の順序
	bool operator<(const Score &other) const {
		return score < other.score;
	}
};

// 従来の編集距離
Score edit_distance(const vector<int> &s, const vector<int> &t, double mut, double indel, int st = 0, int ed = -1){
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

// パラメータ構造体
class Params{
private:
	double mut, indel, dup;
	vector<vector<Score>> unit_to_unit;
	vector<Score> dup_scores;
	vector<Score> indel_scores;

public:
	Params(double m, double i, double d, const vector<vector<int>> &units){
		mut = m;
		indel = i;
		dup = d;

		// (unit, base) -> (unit, base)
		unit_to_unit.resize(units.size(), vector<Score>(units.size()));
		for(int i = 0; i < units.size(); i++){
			for(int j = i; j < units.size(); j++){
				unit_to_unit[i][j] = edit_distance(units[i], units[j], mut, indel);
				unit_to_unit[j][i] = unit_to_unit[i][j];
			}
		}

		// 各ユニットのdup scoreの配列 (単塩基のdup/contはindelとする。ただし塩基がunitとして追加されている場合を除く。)
		dup_scores.resize(units.size(), Score(0.0));
		for(int i = 0; i < units.size(); i++){
			if(i < 4){
				dup_scores[i].set_indel_2(1, indel);
			}else{
				dup_scores[i].set_dup_2(units[i].size(), dup);
			}
		}

		// 各ユニットのindel scoreの配列
		indel_scores.resize(units.size(), Score(0.0));
		for(int i = 0; i < units.size(); i++){
			indel_scores[i].set_indel_2(units[i].size(), indel);
		}
	}

	Score get_unit_to_unit(int i, int j) const {
		return unit_to_unit[i][j];
	}
	Score get_dup(int a) const {
		return dup_scores[a];
	}
	Score get_indel(int a) const {
		return indel_scores[a];
	}

	double mut_val() const {return mut;}
	double indel_val() const {return indel;}
	double dup_val() const {return dup;}

	// Score mut(int b1, int b2){
	// 	if(b1 == b2){
	// 		return Score(0.0);
	// 	}else{
	// 		return Score(mut);
	// 	}
	// }
};

// Score:: 各変異の追加
void Score::set_indel(int len, const Params &par){
	if(score < 0){cerr << "Error: set_indel for minus scores" << endl;}
	score += par.indel_val() * len;
	indel += len;
}
void Score::set_mut(const Params &par){
	if(score < 0){cerr << "Error: set_mut for minus scores" << endl;}
	score += par.mut_val();
	mut++;
}
void Score::set_dup(int len, const Params &par){
	if(score < 0){cerr << "Error: set_dup for minus scores" << endl;}
	score += par.dup_val() * len;
	if(dup.find(len) == dup.end()){
		dup.emplace(len, 1);
	}else{
		dup[len]++;
	}
}

// debug用: 表示
void print_vec3(vector<vector<vector<Score>>> vec){
	cout << "<------ print vec(3) ------>" << endl;
	for(int a = 0; a < vec.size(); a++){
		if(a < 4){continue;}
		cout << "unit " << a << ": " << endl;
		for(vector<Score> inner_vec : vec[a]){
			for(Score score : inner_vec){cout << setw(2) << setfill(' ') << score.get_score() << "  ";}
			cout << endl;
		}
	}
	cout << endl;
}
void print_vec2(vector<vector<Score>> vec){
	cout << "<------ print vec(2) ------>" << endl;
	for(vector<Score> inner_vec : vec){
		for(Score score : inner_vec){cout << setw(2) << setfill(' ') << score.get_score() << "  ";}
		cout << endl;
	}
	cout << endl;
}

// あとでいちいちコピーが生じないように工夫する
// vector<int> string_cut(vector<int> s, int i, int j){
// 	vector<int> vec(j-i);
// 	for(int idx = i; idx < j; idx++){
// 		vec[idx-i] = s[idx];
// 	}
// 	return vec;
// }

// Stage1の計算(だいたい論文通り)
// sとtは方向が逆なだけなので計算方法は同じ
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
	
	// initialization
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
	// 		S[a][i][i] = params.get_indel(a);
	// 	}
	// }
	
	// DP
	for(int j = 2; j <= n; j++){
		for(int i = j-2; i >= 0; i--){
			// S2の更新
			for(int a = 0; a < l; a++){
				for(int h = i+1; h < j; h++){
					Score tmp = min({
						S[a][i][h] + S_eps[h][j],
						S_eps[i][h] + S[a][h][j],
						params.get_dup(a) + S[a][i][h] + S[a][h][j]
					});

					// 元のアルゴリズムだとうまくいかないのでこの部分を追加
					if(j-i <= units[a].size()){
						tmp = min(tmp, edit_distance(s, units[a], params.mut_val(), params.indel_val(), i, j));
					}

					if(S2[a][i][j].get_score() < 0){
						S2[a][i][j] = tmp;
					}else{
						S2[a][i][j] = min(S2[a][i][j], tmp);
					}
				}
			}
			// Sの更新
			for(int a = 0; a < l; a++){
				for(int b = 0; b < l; b++){
					if(S[a][i][j].get_score() < 0){
						S[a][i][j] = params.get_unit_to_unit(a,b) + S2[b][i][j];
					}else{
						S[a][i][j] = min(S[a][i][j], params.get_unit_to_unit(a,b) + S2[b][i][j]);
					}
				}
			}
			// S_epsの更新
			for(int a = 0; a < l; a++){
				if(S_eps[i][j].get_score() < 0){
					S_eps[i][j] = params.get_indel(a) + S[a][i][j];
				}else{
					S_eps[i][j] = min(S_eps[i][j], params.get_indel(a) + S[a][i][j]);
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
	
	vector<vector<Score>> ED(n+1, vector<Score>(m+1, Score(-1.0)));
	vector<vector<vector<Score>>> EDT(l, vector<vector<Score>>(n+1, vector<Score>(m+1, Score(-1.0))));
	
	// S_eps[i,j]: Cost(s[i,j), eps)
	vector<vector<Score>> S_eps(n+1, vector<Score>(n+1, Score(-1.0)));
	vector<vector<vector<Score>>> S(l, vector<vector<Score>>(n+1, vector<Score>(n+1, Score(-1.0))));
	vector<vector<vector<Score>>> S2(l, vector<vector<Score>>(n+1, vector<Score>(n+1, Score(-1.0))));
	vector<vector<Score>> T_eps(m+1, vector<Score>(m+1, Score(-1.0)));
	vector<vector<vector<Score>>> T(l, vector<vector<Score>>(m+1, vector<Score>(m+1, Score(-1.0))));
	vector<vector<vector<Score>>> T2(l, vector<vector<Score>>(m+1, vector<Score>(m+1, Score(-1.0))));
	
	// Stage1
	string_to_unit(s, units, S_eps, S, S2, params);
	// cout << "S_eps "; print_vec2(S_eps);
	// cout << "S "; print_vec3(S);
	// cout << "S2 "; print_vec3(S2);
	string_to_unit(t, units, T_eps, T, T2, params);
	// cout << "T_eps "; print_vec2(T_eps);
	// cout << "T "; print_vec3(T);
	// cout << "T2 "; print_vec3(T2);

	// Stage2
	// Initialization
	ED[0][0].set_zero();
	for(int i = 1; i <= m; i++){
		ED[0][i] = T_eps[0][i];
		ED[1][i] = T[s[0]][0][i];
	}
	for(int i = 1; i <= n; i++){
		ED[i][0] = S_eps[0][i];
		ED[i][1] = S[t[0]][0][i];
	}

	// EDTの初期化(i=1)
	for(int a = 0; a < l; a++){
		for(int j = 2; j <= m; j++){
			for(int h = 1; h < j; h++){
				if(EDT[a][1][j].get_score() < 0){EDT[a][1][j] = ED[1][h] + T[a][h][j];}
				else{EDT[a][1][j] = min(EDT[a][1][j], ED[1][h] + T[a][h][j]);}
			}
		}
	}

	// DP
	for(int j = 2; j <= m; j++){
		for(int i = 2; i <= n; i++){
			// EDTの更新
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
			// EDの更新
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
