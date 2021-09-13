/* Computes the following for every position of the string s:
 *
 * best[0][i] - length of the longest palindrome that ends in s[i]
 * best[1][i] - length of the longest palindrome that starts in s[i]
 * d1[i] - number of odd length palindromes with center in s[i]
 * d2[i] - number of even length palindromes with center in s[i]
 *         (we consider the center of an even length palindrome as 
 *         the rightmost of the two characters in the center)
 *
 * Complexity: O(n)
 */

int best[2][N]; 
int d1[N], d2[N];

void manacher(string &s){
	n = sz(s);
	for(int i = 0; i < n; i++)
		best[0][i] = best[1][i] = 1;
	for (int i = 0, l = 0, r = -1; i < n; i++) {
		int k = (i > r) ? 1 : min(d1[l + r - i], r - i + 1);
		while (0 <= i - k && i + k < n && s[i - k] == s[i + k]) {
			best[1][i-k] = max(best[1][i-k], k << 1 | 1);
			best[0][i+k] = max(best[0][i+k], k << 1 | 1);
			k++;
		}
		d1[i] = k--;
		if (i + k > r) {
			l = i - k;
			r = i + k;
		}
	}
	for (int i = 0, l = 0, r = -1; i < n; i++) {
		int k = (i > r) ? 0 : min(d2[l + r - i + 1], r - i + 1);
		while (0 <= i - k - 1 && i + k < n && s[i - k - 1] == s[i + k]) {
			best[1][i-k-1] = max(best[1][i-k-1], k*2 + 2);
			best[0][i+k] = max(best[0][i+k], k*2 + 2);
			k++;
		}
		d2[i] = k--;
		if (i + k > r) {
			l = i - k - 1;
			r = i + k ;
		}
	}
}