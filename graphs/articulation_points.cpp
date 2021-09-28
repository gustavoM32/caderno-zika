/* Finds the articulation points of an undirected graph. The edges might be parallel or 
 * self-loops without problems.
 *
 * Usage: first add the edges with `addEdge(u, v)`, and then call `findAps(n)`.
 *
 * n - (input) number of vertices.
 * ap - ap[i] is true if the vertex i is an articulation point.
 *
 * Complexity: O(V + E) time and space
 */

int low[N], disc[N];
bool ap[N], t_in;
vector<int> adj[N];

void compute(int v, int root, int p=-1) {
	low[v] = disc[v] = ++t_in;
	int childs = 0;
	for(int u : adj[v]) {
		if(!disc[u]) {
			childs++;
			compute(u, root, v);
			low[v] = min(low[v], low[u]);
			if(low[u] >= disc[v] && (v != root || childs > 1))
				ap[u] = 1;
		} else if(u != p) {
			low[v] = min(low[v], disc[u]);
		}
	}
}

void addEdge(int u, int v) {
	adj[u].pb(v);	
	adj[v].pb(u);	
}

void findAps(int n) {
	for(int i = 0; i < n; i++) {
		low[i] = disc[i] = ap[i] = 0;
	}
	t_in = 0;
	for(int i = 0; i < n; i++) {
        if(!disc[i]) compute(i, i);
	}
}
