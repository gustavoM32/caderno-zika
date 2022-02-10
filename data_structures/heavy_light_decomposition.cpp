/* Data structure to answer queries on paths of a tree. It divides the tree in chains
 * and for each chain saves a data structure that answers the query as if it were done
 * in an array. In this template, this data structure is assumed to be a segment tree
 * but that's not mandatory.
 *
 * Complexity: O(T(n) * logn) time per query/update, where T(n) is the time of
 *              query/update of the inherent data structure.
 */

struct hld {
    int gid, r;
    vector<int> tam, id, p, d, rt;
    vector<vector<int> > adj;
    segTree st;
    hld(int n, int root=0)
       : r(root), gid(0), tam(n), id(n), p(n), d(n), rt(n), adj(n), st(n) {}
    void addEdge(int u, int v) {
        adj[u].pb(v);
        adj[v].pb(u);
    }
    int prec(int v, int par=-1, int depth=0){
        tam[v] = 1;
        p[v] = par;
        d[v] = depth;
        for(int u: adj[v]) if(u != par)
            tam[v] += prec(u, v, depth + 1);
        return tam[v];
    }
    void build(int v, int root){
        id[v] = gid++;
        rt[v] = root;
        if(sz(adj[v]) > 1 && adj[v][0] == p[v])
            swap(adj[v][0], adj[v][1]);
        for(auto &u: adj[v]) if(u != p[v] && tam[u] > tam[adj[v][0]])
            swap(adj[v][0], u);
        for(auto u: adj[v]) if(u != p[v])
            build(u, u == adj[v][0] ? root : u);
    }
    void init(){
        prec(r);
        build(r, r);
    }
    void updateVertex(int u, int x){
        st.update(id[u], x);
    }
    void updateEdge(int u, int v, int x) {
        st.update(max(id[u], id[v]), x);
    }
    // for queries on the edges of the path set for_edge to true
    // this code assumes that the segment tree queries are right-exclusive
    int query(int u, int v, bool for_edge=false){
        auto oper = (ll a, ll b) { // you probably will only need to change this
            return max(a, b);
        };
        ll ans = -INF;
        while(rt[u] != rt[v]){
            if(d[rt[u]] > d[rt[v]]) swap(u, v);
            ans = oper(ans, st.query(id[rt[v]], id[v] + 1));
            v = p[rt[v]];
        }
        int a = id[u], b = id[v];
        ans = oper(ans, st.query(min(a, b) + for_edge, max(a, b) + 1));
        return ans;
    }    
};
