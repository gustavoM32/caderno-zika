#define PROBLEM "https://judge.yosupo.jp/problem/lca"
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
#include "lca.cpp"

int main() {
    fastio;
    int n, q;
    cin >> n >> q;

    for (int i = 1; i < n; i++) {
        int p; cin >> p;
        addEdge(p, i);
    }

    eulerTour(0, 0);

    for (int i = 0; i < q; i++) {
        int u, v;
        cin >> u >> v;
        cout << lca(u, v) << endl;
    }

    return 0;
}