/* 2d iterative segment tree, with point update and range queries.
 * 
 * The operator needs to be commutative for this implementation.
 * 
 * From: https://github.com/mhunicken/icpc-team-notebook-el-vasito/blob/master/data_structures/segment_tree_2d.cpp
 *
 * Time complexity: O(n*m) for building and O(log n * log m) for updates and queries.
 * Space complexity: O(n*m)
 */
struct segTree {
    int n, m;
    vector<vector<ll>> st;
    const ll NEUT = 0; // TODO define neutral element

    // combine two elements, needs to be commutative
    inline ll combine(ll a, ll b) {
        return a + b; // TODO define merge operator
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
        st.assign(2*n, vector<ll>(2*m, NEUT));
    }

    // set position (x, y) to k
    void update(int x, int y, ll k) {
        st[n + x][m + y] = k; // TODO change update operation
        
        for (int j = m + y; j > 1; j >>= 1)
            st[n + x][j >> 1] = combine(st[n + x][j], st[n + x][j ^ 1]);

        for (int i = n + x; i > 1; i >>= 1) for (int j = m + y; j >= 1; j >>= 1)
            st[i >> 1][j] = combine(st[i][j], st[i ^ 1][j]);
    }

    // query in the rectangle (is, js) (ie, je), INCLUSIVE !!!
    ll query(int is, int js, int ie, int je) {
        ll res = NEUT;

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
