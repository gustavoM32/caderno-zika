/* Computes the kmp automata. The i-th state of this automata corresponds
 * to any string which its longest suffix that is a prefix of string s is
 * s[0, i]. The entry kmp[i][j] of the matrix corresponds to the resulting
 * state of transitioning from state i with character START + j
 *
 * START - first character of the alphabet 
 * ALPH - alphabet size
 *
 * Complexity: O(n * ALPH)
 */

const int ALPH = 26;
const char START = 'a';

struct kmpAutomata {
	int pi[N], kmp[N][ALPH];	
	string s;

	int go(int i, int j){
		if(kmp[i][j] != -1) return kmp[i][j];
		int ans;
		if(s[i] == j + START) ans = i + 1;
		else if(i == 0) ans = 0;
		else ans = go(pi[i-1], j);
		return kmp[i][j] = ans;
	}

	kmpAutomata(string _s) s(_s) {
		int n = sz(s);
		int k = 0;
		for(int i = 1; i < n; i++){
			while(k && s[i] != s[k])
				k = pi[k-1];
			if(s[i] == s[k]) k++;
			pi[i] = k;
		}
		for(int i = 0; i < n; i++)
			for(int j = 0; j < ALPH; j++)
				kmp[i][j] = -1;
		for(int i = 0; i < n; i++)
			for(int j = 0; j < ALPH; j++)
				go(i, j);
	}
};