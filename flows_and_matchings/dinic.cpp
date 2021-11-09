/* Fast max flow algorithm.
 *
 * Constructor:
 * dinic(n, s, t)
 * n - number of nodes in the flow network.
 * s - source of the flow network.
 * t - sink of the flow network.
 *
 * Methods:
 * - addEdge(u, v, cap)
 *   adds a directed edge from `u` to `v` with capacity `cap`.
 * - getFlow()
 *   returns the maximum flow of the network.
 *
 * Complexity: In general, the time complexity of getFlow is O(E V^2), but there are better 
 *             upper bounds for bipartite graphs (O(E sqrt(V))) and networks with unit
 *             capacities (O(E sqrt(E))).
 */

const int INF = 1 << 29;

struct dinic {
    ll n, s, t;
    vector<ll> dist, q, work;
    struct edge {
        ll to, rev, f, cap;
    };
    vector<vector<edge> > g;
    dinic(int n, int s, int t) : n(n), s(s), t(t), g(n), dist(n), q(n), work(n) {}
    void addEdge(int u, int v, ll cap) {
        g[u].pb((edge){v, sz(g[v]), 0, cap});
        g[v].pb((edge){u, sz(g[u]) - 1, 0, 0});
    }
    bool bfs() {
        for(int i = 0; i < n; i++) dist[i] = -1;
        dist[s] = 0;
        int qt = 0;
        q[qt++] = s;
        for(int qh = 0; qh < qt; qh++) {
            int u = q[qh];
            for(int i = 0; i < sz(g[u]); i++) {
                edge &e = g[u][i];
                int v = g[u][i].to;
                if(dist[v] < 0 && e.f < e.cap) {
                    dist[v] = dist[u] + 1;
                    q[qt++]=v;
                }
            }
        }
        return dist[t] >= 0;
    }
    ll dfs(int u, ll f) {
        if(u == t) return f;
        for(ll &i = work[u]; i < sz(g[u]); i++) {
            edge &e = g[u][i];
            if(e.cap <= e.f) continue;
            int v = e.to;
            if(dist[v] == dist[u] + 1) {
                ll df = dfs(v, min(f, e.cap - e.f));
                if(df > 0){
                    e.f += df;
                    g[v][e.rev].f -= df;
                    return df;
                }
            }
        }
        return 0;
    }
    ll getFlow() {
        ll res = 0;
        while(bfs()) {
            for(int i = 0; i < n; i++) work[i] = 0;
            while(ll delta = dfs(s, INF))
                res += delta;
        }
        return res;
    }
};