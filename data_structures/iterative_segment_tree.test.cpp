#include "../base_template.cpp"
#include "iterative_segment_tree.cpp"
#define PROBLEM "https://judge.yosupo.jp/problem/point_add_range_sum"

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