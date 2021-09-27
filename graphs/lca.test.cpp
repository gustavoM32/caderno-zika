#include "../base_template.cpp"
#include "lca.cpp"
#define PROBLEM "https://judge.yosupo.jp/problem/lca"

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