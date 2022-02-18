#define PROBLEM "https://judge.yosupo.jp/problem/static_range_sum"
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

#include "mos.cpp"

int main() {
    int n, q;
    cin >> n >> q;
    vector<Query> qs(q);
    vector<ll> v(n);
    for (int i = 0; i < n; i++) {
        cin >> v[i];
    }
    for (int i = 0; i < q; i++) {
        cin >> qs[i].l >> qs[i].r;
        qs[i].r--;
        qs[i].idx = i;
    }
    Mos mos(v);
    vector<ll> ans;
    mos.exec(qs, ans);
    for (ll i : ans) cout << i << "\n";

    return 0;
}