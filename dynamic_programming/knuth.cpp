/* A dp of the form
 *     dp[l][r] = min_{l < m < r}(dp[l][m] + dp[m][r])  + cost(l, r)
 * can be solved in O(n^2) with Knuth opmitization if we have that
 *     opt[l][r - 1] <= opt[l][r] <= opt[l + 1][r]
 * where
 *     dp[l][r] = dp[l][opt[l][r]] + dp[opt[l][r]][r] + cost(l, r).
 *
 * Other sufficient condition (that implies the previous one) is
 *    given a <= b <= c <= d, we have:
 *        - quadrangle inequality: cost(a, c) + cost(b, d) <= cost(a, d) + cost(b, c)
 *        - monotonicity: cost(b, c) <= cost(a, d)
 *
 * Complexity: O(n^2) time
 *             O(n^2) memory
 */

ll knuth(int n) {
    vector<vector<ll>> dp = vector<vector<ll>>(n, vector<ll>(n));
    vector<vector<int>> opt = vector<vector<int>>(n, vector<int>(n));
    for(int k = 0; k <= n; k++) {
        for(int l = 0; l + k <= n; l++) {
            int r = l + k;
            if(k < 2) {
                dp[l][r] = 0; // base case
                opt[l][r] = l;
            }
            dp[l][r] = INF;
            for(int m = opt[l][r - 1]; m <= opt[l + 1][r]; m++) {
                ll cur = dp[l][m] + dp[m][r] + cost(l, r); // must define O(1) cost function
                if(cur < dp[l][r]) {
                    dp[l][r] = cur;
                    opt[l][r] = m;
                }
            }
        }
    }
    return dp[0][n];
}