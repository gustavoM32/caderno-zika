/* Returns the matrix exponentiation
 *
 * base - d x d matrix
 * n - non-negative integer for the exponent
 * res - output param
 *
 * Complexity: O(d^2 * log(n))
 */
void matProd(vector<vector<ll>> &a, vector<vector<ll>> &b, vector<vector<ll>> &c) {
    int n = a.size();
    int m = a[0].size();
    int p = b[0].size();
    vector<vector<ll>> res(n, vector<ll>(p, 0));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < p; j++) {
            for (int k = 0; k < m; k++) {
                res[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    c = res;
}

void matExp(vector<vector<ll>> &base, ll n, vector<vector<ll>> &res) {
    int d = base.size();
    for (int i = 0; i < d; i++) {
        for (int j = 0; j < d; j++) {
            res[i][j] = ll(i == j);
        }
    }
    while (n > 0) {
        if (n % 2 == 1) matProd(res, base, res);
        matProd(base, base, base);
        n /= 2;
    }
}
