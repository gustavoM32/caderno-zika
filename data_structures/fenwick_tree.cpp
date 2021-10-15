/* Fenwick tree, with point update and range queries.
 * 
 * The operator needs to be commutative for this implementation.
 *
 * Time complexity: O(n log(n)) for building and O(log n) for updates and queries.
 * Space complexity: O(n)
 */
struct Fenwick {
    vector<ll> bit;
    int n;

    void build(vector<ll> &v) {
        for (size_t i = 0; i < v.size(); i++) {
            update(i, v[i]);
        }
    }

    Fenwick() {}

    Fenwick(int n) : n(n) {
        resize(n);
    }

    Fenwick(vector<ll> &v) : Fenwick(v.size()) {
        build(v);
    }

    void resize(int n) {
        bit.assign(n, 0);
    }

    ll query(int r) {
        ll ret = 0;
        for (; r >= 0; r = (r & (r + 1)) - 1) {
            ret += bit[r];
        }
        return ret;
    }

    ll query(int l, int r) {
        return query(r) - query(l-1);
    }

    void update(int i, ll v) {
        for (; i < n; i = i | (i + 1)) {
            bit[i] += v;
        }
    }
};
