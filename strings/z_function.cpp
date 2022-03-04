/* Computes z function, where z[i] is the length of the longest prefix of s[i, n)
 * that is also a prefix of s[0, n). For convention, this template assumes z[0] = 0
 *
 * Complexity: O(n)
 */

vector<int> z_function(string s) {
	int n = sz(s);
	vector<int> z(n);
	for(int i = 1, l = 0, r = 0; i < n; i++) {
		if(i <= r)
			z[i] = min(r - i + 1, z[i - l]);
		while(i + z[i] < n && s[z[i]] == s[i + z[i]])
			z[i]++;
		if(i + z[i] - 1 > r) {
			l = i;
			r = i + z[i] - 1;
		}
	}
	return z;
}