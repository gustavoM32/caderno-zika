/* Computes the Strongly Connected Components of a graph
 *
 * n - (input) number of nodes of the graph
 * adj - (input) vector of adjacency of the graph
 * cont - number of SCCs of the graph
 * SCC[i] - list of nodes inside the i-th SCC
 * belong[i] - the index of the SCC in which the node i is in
 * adjSCC - vector of adjacency of the compressed graph
 *
 * Complexity: O(n) to compute cont, SCC and belong
 *             O(n log(n)) to compute adjSCC 
 */

int n, temp, cont;
vector<int> adj[N], SCC[N], adjSCC[N];
bool vis[N];
int low[N], disc[N], belong[N];
stack<int> S;
 
void dfs(int v) {
	disc[v] = low[v] = ++temp;
	S.push(v);
	vis[v] = 1;
	for(int u: adj[v]) {
		if(!disc[u])
			dfs(u);
		if(vis[u])
			low[v] = min(low[u], low[v]);
	}
	if(disc[v] == low[v]) {
		while(1) {
			int u = S.top();
			S.pop();
			vis[u] = 0;
			belong[u] = cont;
			SCC[cont].pb(u);
			if(u == v) break;
		}
		cont++;
	}
}
 
void tarjan() {
	for(int i = 0; i < n; i++)
		if(!disc[i])
			dfs(i);
	for(int i = 0; i < cont; i++) {
		for(int u: SCC[i])
			for(int v: adj[u]) if(belong[v] != i) {
				adjSCC[i].pb(belong[v]);
			}
		sort(adjSCC[i].begin(), adjSCC[i].end());
		adjSCC[i].resize(unique(adjSCC[i].begin(), adjSCC[i].end()) - adjSCC[i].begin());
	}
}