/* Modular inverse for numbers 1 to n mod m.
 *
 * Complexity: O(n)
 */
void inverse_1n(int n, int m) {
    vector<int> inv(n+1);
    inv[1] = 1;
    for (int i = 2; i <= n; i++) {
        inv[i] = m - ((m/i) * inv[m%i]) % m;
    }

}