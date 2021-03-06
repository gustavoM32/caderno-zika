/* Simple iterative segment tree, with point update and range queries.
 * 
 * The operator needs to be commutative for this implementation.
 *
 * Time complexity: O(n) for building and O(log n) for updates and queries.
 * Space complexity: O(n)
 */
struct segTree {
    int n;
    vector<ll> st;
    const ll NEUT = 0; // TODO define neutral element

    // combine two elements, doesn't need to be commutative
    inline ll combine(ll a, ll b) {
        return a + b; // TODO define merge operator
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

    segTree() {}

    // initialize with neutral elements
    segTree(int n) {
        resize(n);
    }

    // initialize with vector
    segTree(vector<ll> &v) : segTree(v.size()) {
        build(v);
    }

    void resize(int s) {
        n = s;
        st.assign(2*s, NEUT);
    }

    // add x to position i
    void update(int i, ll x) {
        st[i += n] += x; // TODO change update operation
        while (i > 1) {
            i >>= 1;
            st[i] = combine(st[i << 1], st[i << 1 | 1]);
        }
    }

    // query from l to r, inclusive
    ll query(int l, int r) {
        ll resl = NEUT, resr = NEUT;
        for (l += n, r += n+1; l < r; l >>= 1, r >>= 1) {
            if (l & 1) resl = combine(resl, st[l++]);
            if (r & 1) resr = combine(st[--r], resr);
        }

        return combine(resl, resr);
    }
};
