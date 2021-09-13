/* Random number generator
 *
 * Complexity: O(1)
 */
mt19937 rng(chrono::steady_clock::now().time_since_epoch().count());
int x = rng() % n; // to generante random number
