/* Finds the bridges of an undirected graph. The edges might be parallel or self-loops
 * without problems.
 *
 * Usage: first add the edges with `addEdge(u, v)`, and then call `findBridges(n)`
 *
 * n - (input) number of vertices.
 * bridge - bridge[i] is true if the i-th edge is a bridge.
 *
 * Complexity: O(V + E) time and space
 */

vector<pair<int, int>> edge;
int low[N], disc[N];
bool bridge[N], t_in;
vector<int> adj[N];

void compute(int v, int p=-1) {
	low[v] = disc[v] = ++t_in;
	for(int e : adj[v]) {
		int u = edge[e].first ^ edge[e].second ^ v;
		if(!disc[u]) {
			childs++;
			compute(u, e);
			low[v] = min(low[v], low[u]);
			if(low[u] > disc[v]) bridge[e] = 1;
		} else if(e != p) {
			low[v] = min(low[v], disc[u]);
		}
	}
}

void addEdge(int u, int v) {
	adj[u].pb(sz(edge));	
	adj[v].pb(sz(edge));	
	edge.pb({u, v});
}

void findBridges(int n) {
	for(int i = 0; i < n; i++) {
		low[i] = disc[i] = 0;
	}
	for(int i = 0; i < sz(edges); i++) {
		bridges[i] = 0;
	}
	t_in = 0;
	for(int i = 0; i < n; i++) {
        if(!disc[i]) compute(i);
	}
}
