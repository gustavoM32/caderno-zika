#define PROBLEM "https://judge.yosupo.jp/problem/point_add_range_sum"
/**********/
#include <bits/stdc++.h>
using namespace std;

#define fastio ios_base::sync_with_stdio(0);cin.tie(0)
#define pb push_back
#define mp make_pair
#define sz(x) int(x.size())
#define trace(x) cerr << #x << ": " << x <<endl;

typedef long long ll;

const ll N = 1e6;
const ll INF = 1LL << 61;
const ll MOD = 1e9 + 7;
/**********/
#include "iterative_segment_tree.cpp"

int main() {
    int n, q;
    cin >> n >> q;

    vector<ll> v(n);

    for (int i = 0; i < n; i++) {
        cin >> v[i];
    }

    segTree st(v);
    
    for (int i = 0; i < q; i++) {
        int op;
        cin >> op;
        if (op == 0) {
            int p, x;
            cin >> p >> x;
            st.update(p, x);
        } else {
            int l, r;
            cin >> l >> r;
            cout << st.query(l, r-1) << "\n";
        }
    }

    return 0;
}