/* Returns the modular multiplicative inverse of a number or -1 if gcd(a, m) != 1
 *
 * a, m - positive integers
 *
 * Complexity: O(log(m))
 */
ll gcdExt(ll a, ll b, ll &x, ll &y) {
	if (b == 0) {
		x = 1;
		y = 0;
		return a;
	}
	ll x1, y1;
	ll res = gcdExt(b, a % b, x1, y1);
	x = y1;
	y = x1 - y1 * (a / b);
	return res;
}

ll inv(ll a, ll m = MOD) {
	ll x, y;
	ll g = gcdExt(a, m, x, y);
	if (g != 1) return -1;
	return (m + x % m) % m;
}