/* Solves the bipartite matching problem for color classes X and Y.
 *
 * n - (input) size of color class X
 * m - (input) size of color class Y
 * g - (input) g[x] contains the neighbours of x \in X 
 * mat - mat[y] is the vertex matched with y \in Y or -1 if there's no such vertex
 * arcs - set of arcs of the bipartite matching. Each arc is a pair such that the first
          element is in X and the second in Y.  
 *
 * Complexity: O(V * (V + E)) where V is the number of vertex and E the number of edges
                              in the bipartite graph  
 */

vector<int> g[N];
int mat[N];
bool vis[N];
int n, m;
 
int match(int x) {
	if(vis[x]) return 0;
	vis[x] = true;
	for(int y: g[x]) if(mat[y] < 0 || match(mat[y])) {
		mat[y] = x;
		return 1;
	}
	return 0;
}
 
vector<pair<int, int> > maxMatching() {
	vector<pair<int, int> > arcs;
	for(int i = 0; i < m; i++)
		mat[i] = -1;
	for(int i = 0; i < n; i++) {
		for(int j = 0; j < n; j++)
			vis[j] = 0;
		match(i);
	}
	for(int i = 0; i < m; i++)
		if(mat[i] >= 0) arcs.pb({mat[i], i});
	return arcs;
}