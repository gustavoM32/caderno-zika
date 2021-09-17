/* Extended Euclidean algorithm. Returns gcd(a, b) and set the parameters
 * x and y to numbers such that ax + by = gcd(a, b).
 *
 * Time complexity: O(log(min(a, b)))
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
