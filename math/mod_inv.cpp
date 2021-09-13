/* Returns the modular multiplicative inverse of a number
 *
 * a - integer
 * m - prime module
 *
 * Complexity: O(log(m))
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

ll inv(ll a, ll m = MOD) {
	return modexp(a, m-2, m);
}