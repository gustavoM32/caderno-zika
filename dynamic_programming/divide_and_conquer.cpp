/* A dp of the form
 *     dp[i][j] = min_{k < j}(dp[i - 1][k]  + cost(k, j))
 * can be solved in O(m n logn) with divide and conquer optimization if we have that
 *     opt[i][j] <= opt[i][j + 1]
 * where
 *     dp[i][j] = dp[i - 1][opt[i][j]]  + cost(opt[i][j], j).
 *
 * Complexity: O(m n logn) time (for a partition in m subarrays of an array of size n)
 *             O(n) memory
 */

ll dp[N][2];

void go(int k, int l, int r, int optl, int optr) {
    if(l > r) return;
    int opt, m = (l + r) >> 1;
    dp[m][k&1] = INF;
    for(int i = optl; i <= min(m, optr + 1); i++) {
        ll cur = dp[i][~k&1] + cost(i, m); // must define O(1) cost function
        if(cur < dp[m][k&1]) {
            dp[m][k&1] = cur;
            opt = i;
        }
    }
    go(k, l, m - 1, optl, opt);
    go(k, m + 1, r, opt, optr);
}

ll dc(int n, int m) {
    dp[0][0] = dp[0][1] = 0;
    for(int i = 0; i < n; i++) dp[i][0] = INF;
    for(int i = 1; i <= m; i++) go(i, i , n, i - 1, n - 1);
    return dp[n - 1][m & 1];
}