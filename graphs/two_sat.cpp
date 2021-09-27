/* Returns true if a given set of conditions of the form a v b is
 * satisfiable. NOTE: you must clear adj and adjt between testcases.
 *
 * n - (input) number of atomic propositions
 * assign - one valid assignation of truth values
 *
 * Complexity: O(n)
 */

// Useful when negation of atomic expression is given as a negative number
int conv(int x) {
	if(x > 0) return (x - 1) << 1;
	else return (-x - 1) << 1 | 1;
}
 
vector<int> adj[N], adjt[N];
 
void addEdge(int u, int v) {
	adj[u].pb(v);
	adjt[v].pb(u);
}
 
// Add the condition a v b
void addProp(int a, int b) {
	a = conv(a), b = conv(b);
	addEdge(a^1, b);
	addEdge(b^1, a);
}
 
bool vis[N];
int comp[N];
stack<int> order;
 
void toposort(int v) {
	vis[v] = 1;
	for(int u: adj[v]) if(!vis[u])
		toposort(u);
	order.push(v);
}
 
void mark(int v, int id) {
	comp[v] = id;
	for(int u: adjt[v]) if(comp[u] == -1)
		mark(u, id);
}

bool assign[N];

// Call this function after adding al conditions with addProp
bool twoSat(int n) {
	for(int i = 0; i < (n << 1); i++) {
		comp[i] = -1;
		vis[i] = 0;
	}
	for(int i = 0; i < (n << 1); i++) if(!vis[i])
		toposort(i);
	int cont = 0;
	for(int i = 0; i < (n << 1); i++) {
		int v = order.top();
		order.pop();
		if(comp[v] == -1)
			mark(v, cont++);
	}
	for(int i = 0; i < n; i++) {
		if(comp[i << 1] == comp[i << 1 | 1])
			return 0;
		assign[i] = comp[i << 1] > comp[i << 1 | 1];
	}
	return 1;
}