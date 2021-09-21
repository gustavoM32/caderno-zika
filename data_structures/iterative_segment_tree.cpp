/* Simple iterative segment tree, with point update and range queries.
 * 
 * The operator needs the be commutative for this implementation.
 *
 * Time complexity: O(n) for building and O(log n) for updates and queries.
 * Space complexity: O(n)
 */
struct segTree {
    int n;
    vector<ll> st;

    // combine two elements
    inline ll combine(ll a, ll b) {
        return a + b;
    }

    // build the tree with vector v
    void build(vector<ll> &v) {
        for (int i = 0; i < n; i++) {
            st[n + i] = v[i];
        }

        for (int i = n-1; i >= 1; i--) {
            st[i] = combine(st[i << 1], st[i << 1 | 1]);
        }
    }

    public:

    // initialize with zeroes
    segTree(int n) : n(n) {
        st.assign(2*n, 0);
    }

    // initialize with vector
    segTree(vector<ll> &v) : segTree(v.size()) {
        build(v);
    }

    // update position i with value x
    void update(int i, ll x) {
        i += n;
        st[i] = x;
        while (i > 1) {
            st[i >> 1] = combine(st[i], st[i ^ 1]);
            i >>= 1;
        }
    }

    // query from l to r, inclusive
    ll query(int l, int r) {
        ll res = 0;
        l += n;
        r += n+1;
        while (l < r) {
            if (l & 1) res = combine(res, st[l++]);
            if (r & 1) res = combine(res, st[--r]);

            l >>= 1;
            r >>= 1;
        }

        return res;
    }
};