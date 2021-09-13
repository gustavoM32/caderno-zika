/* Computes pi, where pi[i] is the length of the longest (proper) prefix of s[0, i]
 * that is also a (proper) suffix of s[0, i]
 *
 * Complexity: O(n)
 */

int pi[N];

void kmp(string& s) {
	int n = sz(s), k = 0;
	for(int i = 1; i < n; i++){
		while(k && s[i] != s[k])
            k = pi[k-1];
		if(s[i] == s[k]) k++;
		pi[i] = k;
	}
}