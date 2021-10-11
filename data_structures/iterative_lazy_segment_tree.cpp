/* Iterative segment tree with lazy propagation. Supports range updates and queries.
 * This example is for querying the maximum value in a range and updating with sum in
 * a range.
 *
 * Changes must be done in the struct node operators nad NEUT value. For more complicated
 * lazy values, change the `apply` method in segTree too.
 *
 * Time complexity: O(n) for building and O(log n) for updates and queries.
 * Space complexity: O(n)
 */

const ll NEUT = -INF;

struct node {
    ll val;
    node() : val(0) {} // initial value
    node(ll val) : val(val) {}
    // combine two nodes
    node operator+(const node& other) {
        return node(max(val, other.val));
    }
    // update a node by the lazy value
    void operator+=(ll x) {
        val += x;
    }
};

struct segTree {
    int n, h;
    vector<ll> d;
    vector<node> t;
    segTree(int n) : n(n), t(n << 1) {
        d.resize(n, 0);
        h = sizeof(int) * 8 - __builtin_clz(n);
    }
    void apply(int p, ll x) {
        t[p] += x;
        if(p < n) d[p] += x;
    }
    void push(int p) {
        for(int s = h; s > 0; s--) {
            int i = p >> s;
            if(d[i]) {
                apply(i << 1, d[i]);
                apply(i << 1 | 1, d[i]);
                d[i] = 0;
            }
        }
    }
    void build(int p) {
        for(; p >>= 1;) {
            t[p] = t[p << 1] + t[p << 1 | 1];
            t[p] += d[p];
        }
    }
    void update(int l, int r, ll x){
        l += n, r += n;
        int l0 = l, r0 = r;
        push(l);
        push(r - 1);
        for(; l < r; l >>= 1, r >>= 1) {
            if(l & 1) apply(l++, x);
            if(r & 1) apply(--r, x);
        }
        build(l0);
        build(r0 - 1);
    }
    node query(int l, int r) {
        l += n, r += n;
        push(l);
        push(r - 1);
        node ans(NEUT);
        for(; l < r; l >>= 1, r >>= 1){
            if(l & 1) ans = ans + t[l++];
            if(r & 1) ans = ans + t[--r];
        }
        return ans;
    }
};