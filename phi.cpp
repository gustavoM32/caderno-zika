/* Count how many numbers 1 to n are relatively prime with n
 *
 * n - positive integer
 *
 * Complexity: O(sqrt(n))
 */
ll phi(ll n) {
    ll res = n;
    for (ll p = 2; p * p; p++) {
        if (n % p == 0) {
            while (n % p == 0) n /= p;
            res -= res / p;
        }
    }
    if (n > 1) res -= res / n;
    return res;
}
