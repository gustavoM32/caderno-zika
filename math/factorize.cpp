/* Gets prime factors of a number
 *
 * a - positve integer
 * f - output vector for the factors
 *
 * Complexity: O(sqrt(n))
 */
void getFactors(ll n, vector<ll> &factors) {
    for (ll p = 2; p * p <= n; p++) {
        while (n % p == 0) {
            factors.push_back(p);
            n /= p;
        }
    }
    if (n > 1) factors.push_back(n);
}

