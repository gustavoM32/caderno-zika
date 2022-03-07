/* Sparse table. Useful for queries of idempotent functions in a range.
 * Change the oper method accordingly.
 *
 * Time complexity: 
 *   - O(n logn) for building the structure
 *   - O(1) for queries of idempotent functions.
 * Space complexity: O(n logn)
 */

struct sparseTable {
    int n, logn;
    vector<vector<ll>> t;
    int log_floor(int tam) {
        return 31 - __builtin_clz(tam);
    }
    ll oper(ll a, ll b) { // example with minimum in range
        return min(a, b);
    }
    sparseTable() {}
    sparseTable(vector<ll>& v) {
        n = sz(v);
        logn = log_floor(n) + 1;
        t.resize(logn);
        t[0].resize(n);
        for(int i = 0; i < n; i++)
            t[0][i] = v[i];
        for(int k = 1; k < logn; k++) {
            t[k].resize(n);
            for(int i = 0; i + (1 << k) <= n; i++)
                t[k][i] = oper(t[k - 1][i], t[k - 1][i + (1 << (k - 1))]);
        }
    }
    ll que(int l, int r) { // queries in the semi-open interval [l, r)
        int k = log_floor(r - l);
        return oper(t[k][l], t[k][r - (1 << k)]);
    }
};