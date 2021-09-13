/* Random number generator
 *
 * Complexity: O(1)
 */
mt19937 rng(chrono::steady_clock::now().time_since_epoch().count());
int x = rng() % n; // to generante random number

// to generate a uniform random number in the range [a, b]
ll random(ll a, ll b) { 
    return uniform_int_distribution<ll> (a, b) (rng); 
}