/* Get divisors of a number
 *
 * a - number
 * d - output param to store the divisors (needs to be empty)
 *
 * Complexity: O(n)
 */
void getDivisors(ll a, vector<ll> &d) {
    for (ll i = 1; i * i <= a; i++) {
        if (a % i == 0) {
            d.push_back(i);
            if (i * i != a) d.push_back(a / i);
        }
    }
}
