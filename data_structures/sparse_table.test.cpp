#define PROBLEM "https://judge.yosupo.jp/problem/staticrmq"
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

const ll N = 1e6;
const ll INF = 1LL << 61;
const ll MOD = 1e9 + 7;

#include "sparse_table.cpp"

int main() {
    int n, q;
    cin >> n >> q;

    vector<ll> v(n);

    for (int i = 0; i < n; i++) {
        cin >> v[i];
    }

    sparseTable spt(v);
    
    for (int i = 0; i < q; i++) {
        int l, r; cin >> l >> r;
        l--;
        cout << spt.que(l, r) << endl;
    }

    return 0;
}