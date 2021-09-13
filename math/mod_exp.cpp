/* Logarithmic modular exponentiation
 *
 * b - base
 * e - exponent
 * m - module
 *
 * Complexity: O(log(e))
 */
ll modexp(ll b, ll e, ll m = MOD) {
	ll res = 1;
	b %= m;
	while (e > 0) {
		if (e & 1) res = (res * b) % m;
		b = (b * b) % m;
		e /= 2;
	}
	return res;
}
