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

/*
	Returns true if n is composite

	n - positive integer
	a - fermat base 2 <= a <= n - 2
	d, s - greatest s such that = 2^s * d
*/
bool checkComposite(ll n, ll a, ll d, int s) {
	ll x = modexp(a, d, n);
	if (x == 1 || x == n - 1) return false;
	for (int i = 1; i < s; i++) {
		x = (x * x) % n;
		if (x == n - 1) return false;
	}
	return true;
}

/* Returns true if a number is prime. Deterministic for 64-bit integers.
 *
 * n - positive integer
 *
 * Complexity: O(log(n)) the constant is at least 12
 */
bool millerRabin(ll n) {
	if (n < 2) return false;

	int r = 0;
	ll d = n - 1;
	while ((d & 1) == 0) {
		d >>= 1;
		r++;
	}

	for (int a : {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37}) {
		if (n == a) return true;
		if (checkComposite(n, a, d, r)) return false;
	}
	return true;
} 
