#define PROBLEM "https://codeforces.com/problemset/problem/869/E" // TLE, maybe choose another problem
#define IGNORE // incompatible judge
/**********/
#include <bits/stdc++.h>
using namespace std;

#define fastio ios_base::sync_with_stdio(0);cin.tie(0)
#define pb push_back
#define mp make_pair
#define sz(x) int(x.size())
#define trace(x) cerr << #x << ": " << x <<endl;
/**********/
typedef long long ll;

const ll N = 2505;
const ll INF = 1LL << 61;
const ll MOD = 1e9 + 7;

struct segTree {
    int n, m;
    vector<vector<ll>> st;
    const ll NEUT = -INF; // TODO define neutral element

    // combine two elements, needs to be commutative
    inline ll combine(ll a, ll b) {
        return max(a, b); // TODO define merge operator
    }

    // build the tree with matriz mat
    void build(vector<vector<ll>> &mat) {
        for (int i = 0; i < n; i++) for (int j = 0; j < m; j++)
            st[n + i][m + j] = mat[i][j];

        for (int i = 0; i < n; i++) for (int j = m-1; j >= 1; j--)
            st[n + i][j] = combine(st[n + i][j << 1], st[n + i][j << 1 | 1]);

        for (int i = n-1; i >= 1; i--) for (int j = 0; j < 2*m; j++)
            st[i][j] = combine(st[i << 1][j], st[i << 1 | 1][j]);
    }

    public:

    segTree() {}

    // initialize with neutral elements
    segTree(int n, int m) {
        resize(n, m);
    }

    // initialize with matrix
    segTree(vector<vector<ll>> &m) : segTree(m.size(), m.front().size()) {
        build(m);
    }

    void resize(int new_n, int new_m) {
        n = new_n;
        m = new_m;
        st.assign(2*n, vector<ll>(2*m, NEUT)); // uses NEUT
    }

    // set position (x, y) to k
    void update(int x, int y, ll k) {
        st[n + x][m + y] = k;
        
        for (int j = m + y; j > 1; j >>= 1)
            st[n + x][j >> 1] = combine(st[n + x][j], st[n + x][j ^ 1]);

        for (int i = n + x; i > 1; i >>= 1) for (int j = m + y; j >= 1; j >>= 1)
            st[i >> 1][j] = combine(st[i][j], st[i ^ 1][j]);
    }

    // query in the rectangle (is, js) (ie, je), INCLUSIVE !!!
    ll query(int is, int js, int ie, int je) {
        ll res = NEUT; // uses NEUT

        for (int i0 = n + is, i1 = n + ie + 1; i0 < i1; i0 >>= 1, i1 >>= 1) {
            ll t[2], q = 0;
            if (i0 & 1) t[q++] = i0++;
            if (i1 & 1) t[q++] = --i1;
            for (int k = 0; k < q; k++) {
                for (int j0 = m + js, j1 = m + je + 1; j0 < j1; j0 >>= 1, j1 >>= 1) {
                    if (j0 & 1) res = combine(res, st[t[k]][j0++]);
                    if (j1 & 1) res = combine(res, st[t[k]][--j1]);
                }
            }
        }

        return res;
    }
};
 
void solve() {
    vector<pair<ll, pair<int, int>>> ps;
    vector<pair<ll, pair<pair<int, int>, int>>> qs;
    int n, m, q;
    cin >> n >> q;
    m = n;
    vector<ll> v(n);
    vector<ll> res(q);
    
    segTree st(n, m);
 
    for (int i = 0; i < n; i++) {
        cin >> v[i];
    }
 
    for (int i = 0; i < n; i++) {
        ll sum = 0;
        for (int j = i; j < n; j++) {
            sum += v[j];
            ps.push_back({sum, {i, j}});
        }
    }
 
    for (int i = 0; i < q; i++) {
        int l, r;
        ll u;
        cin >> l >> r >> u;
        qs.push_back({u, {{l-1, r-1}, i}});
    }
 
    sort(ps.begin(), ps.end());
    sort(qs.begin(), qs.end());
 
    int ip = 0;
 
    for (auto q : qs) {
        ll u = q.first;
        int l = q.second.first.first;
        int r = q.second.first.second;
        int id = q.second.second;
        while (ip < ps.size() && ps[ip].first <= u) {
            ll pu = ps[ip].first;
            int i = ps[ip].second.first;
            int j = ps[ip].second.second;
            st.update(i, j, pu);
            ip++;
        }
        res[id] = st.query(l, l, r, r);
    }
 
    for (int i = 0; i < q; i++) {
        if (res[i] == -INF) cout << "NONE\n";
        else cout << res[i] << "\n";
    }
}
 
int main() {
    fastio;
    solve();
    return 0;
}
