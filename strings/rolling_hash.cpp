/* Computes the suffixes of a string s. Can answer queries of the
 * any substring of s.
 *
 * START - first character of the alphabet 
 * K - number of prime bases
 * P - array of prime bases
 * 
 * Complexity: O(n) precomputation, O(1) queries
 */

const ll MOD = 1e9 + 7;
const ll K = 2;
const ll P[K] = {29, 31};
const char START = 'a';

array<ll, K> pot[N];
 
ll mul(ll a, ll b) {
	return a * b % MOD;
}

ll sum(ll a, ll b) {
	a += b;
	if(a >= MOD) a -= MOD;
	return a;
}

// Must be called once in the beginning of the program to compute 
// powers of the prime bases
void prec() {
	for(int k = 0; k < K; k++)
		pot[0][k] = 1;
	for(int i = 1; i < N; i++)
		for(int k = 0; k < K; k++)
			pot[i][k] = mul(pot[i-1][k], P[k]);
}
 
struct hsh {
	array<ll, K> suf[N];
	hsh() {}
	hsh(string &s) {
		for(int i = sz(s)-1; i >= 0; i--)
			for(int k = 0; k < K; k++)
				suf[i][k] = sum(mul(suf[i+1][k], P[k]), s[i] - START + 1);
	}
	// Queries the hashing of the substring s[l, r)
	array<ll, K> que(int l, int r) {
		array<ll, K> cur;
		for(int k = 0; k < K; k++)
			cur[k] = sum(suf[l][k], MOD - mul(suf[r][k], pot[r-l][k]));
		return cur;
	}
};