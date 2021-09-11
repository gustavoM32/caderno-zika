/* Greatest common divisor
 *
 * a, b - non-negative integers
 *
 * Complexity: O(log(min(a, b)))
 */
ll gcd(ll a, ll b) {
	return b ? gcd(b, a % b) : a;
}
