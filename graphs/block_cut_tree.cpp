/* Builds the block-cut tree of a given graph. The block-cut tree is a graph
 * decomposition in biconnected components (blocks) and articulation points (cuts).
 * Every block is uniquely defined by a set of edges. This template considers
 * the the articulation point only belongs to the cut.
 * Usage: add undirected edges with addEdge(u, v) and then call build(n), when n is
 * the total number of vertices in the graph. WARNING: make sure to do not add
 * parallel edges. Also, notice that you have to clear adj, adj_bct and edge_cont
 * by yourself between different test cases.
 *
 * n - (input) number of vertices.
 * belong - belong[i] is the component in which the vertex i is in, either a block
            or a cut.
 * edge_belong - edge_belong[i] is the block in which the edge i is in.
 * ap - ap[i] is true if the vertex i is an articulation point. WARNING: this template
        considers isolated vertices as articulation points.
 * adj_bct - the list of adjacency of the block-cut tree
 *
 * Complexity: O(n) to find all biconnected components and articulation points.
               O(n log n) to build the block-cut tree.
 */

int edge_cont, comp, temp, stp;
pair<int, int> edge[N];
int belong[N], low[N], disc[N], st[N], edge_belong[N];
bool ap[N];
vector<int> adj[N], adj_bct[N];

 
void bcc(int v, int p=-1) {
	low[v] = disc[v] = ++temp;
	int child = 0;
	if(sz(adj[v]) == 0) ap[v] = 1;
	for(int ed: adj[v]) {
		int u = edge[ed].first ^ edge[ed].second ^ v;
		if(!disc[u]) {
			child++;
			st[stp++] = ed;
			bcc(u, v);
			low[v] = min(low[v], low[u]);
			if(low[u] >= disc[v]) {
				if(p != -1 || child > 1) ap[v] = 1;
				while(1) {
					int e = st[--stp];
					edge_belong[e] = comp;
					belong[edge[e].first] = comp;
					belong[edge[e].second] = comp;
					if(e == ed) break;
				}
				comp++;
			}
		} else if(u != p) {
			low[v] = min(low[v], disc[u]);
			if(disc[u] < disc[v])
				st[stp++] = ed;
		}
	}
}
 
 
void addEdge(int u, int v) {
	adj[u].pb(edge_cont);	
	adj[v].pb(edge_cont);	
	edge[edge_cont++] = {u, v};
}
 
void build(int n) {
	for(int i = 0; i < n; i++) {
		low[i] = disc[i] = ap[i] = 0;
	}
	temp = comp = stp = 0;
	for(int i = 0; i < n; i++) 
        if(!disc[i])
		    bcc(i);
	for(int u = 0; u < n; i++) if(ap[u]) {
		for(int e: adj[u]) {
			adj_bct[comp].pb(edge_belong[e]);
			adj_bct[edge_belong[e]].pb(comp);
		}
		belong[u] = comp;
		comp++;
	}
	for(int i = 0; i < comp; i++) {
		sort(all(adj_bct[i]));
		adj_bct[i].resize(unique(all(adj_bct[i])) - adj_bct[i].begin());
	}
}