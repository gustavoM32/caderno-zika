/* Returns the start index of the minimum cyclic string
 * Uses Lyndon Factorization
 *
 * Complexity: O(n)
 */
int min_cyclic_string(string s) {
    s += s;
    int n = s.size();
    int i = 0, res = 0;

    while (i < n / 2) {
        res = i;
        int j = i + 1, k = i;

        while (j < n && s[k] <= s[j]) {
            if (s[k] < s[j]) {
                k = i;
            } else {
                k++;
            }
            j++;
        }

        while (i <= k) {
            i += j - k;
        }
    }

    return res;
}
