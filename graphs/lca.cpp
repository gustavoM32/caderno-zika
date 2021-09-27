/* Binary lifting algorithm to compute the lowest common ancestor of two nodes in a tree
 *
 * To use this template, first you need to add the (undirected) edges using `addEdge` and,
 * after all the edges has been added, call `eulerTour(root, root)` where `root` is the 
 * root of the tree. Then you can use `lca(u, v)` to fin the lowest common ancestor of two
 * nodes `u` and `v`.
 * 
 * LOGN - logarithm in base 2 of N
 * anc[u][j] - (output) 2^j-th ancestor o u
 *
 * Time complexity: O(n log n) for precomputing and O(log n) per query.
 * Space complexity: O(n log n)
 */

const int LOGN = 20;

int anc[N][LOG], tam[N], tin[N];
vi adj[N];
int node_id;

void addEdge(int u, int v) {
    adj[u].pb(v);
    adj[v].pb(u);
}

// Call this with eulerToor(root, root)
int eulerTour(int u, int p) {
    anc[u][0] = p;
    tam[u] = 1;
    tin[u] = node_id++;
    for(int i = 1; i < LOGN; i++)
        anc[u][i] = anc[anc[u][i-1]][i-1];
    for(int v: adj[u]) if(v != p)
        tam[u] += dfs(v, u);
    return tam[u];
}

// Check if u is ancestor of v in O(1)
bool isAnc(int u, int v) {
    return tin[u] <= tin[v] && tin[v] < tin[u] + tam[u];
}

int lca(int u, int v) {
    if(isAnc(u, v)) return u;
    if(isAnc(v, u)) return v;
    for(int i = LOGN - 1; i >= 0; i--)
        if(!isAnc(anc[u][i], v)) u = anc[u][i];
    return anc[u][0];
}