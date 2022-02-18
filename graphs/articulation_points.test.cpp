#define PROBLEM "https://onlinejudge.org/index.php?option=com_onlinejudge&Itemid=8&category=24&page=show_problem&problem=251"
#define IGNORE // incompatible judge
/**********/
#include <bits/stdc++.h>
using namespace std;

#define fastio ios_base::sync_with_stdio(0);cin.tie(0)
#define pb push_back
#define mp make_pair
#define sz(x) int(x.size())
#define trace(x) cerr << #x << ": " << x <<endl;
/**********/
typedef long long ll;

const ll N = 1e6;
const ll INF = 1LL << 61;
const ll MOD = 1e9 + 7;

int low[N], disc[N], t_in;
bool ap[N];
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
				ap[v] = 1;
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

// Must be call between testcases
void clearGraph(int n) {
	for(int i = 0; i < n; i++)
		adj[i].clear();
}

int main() {
	fastio;
	int n;
	while(cin >> n) {
		if(n == 0) break;
		string s;
		while(1) {
			cin >> ws;
			getline(cin, s);
			istringstream iss(s);
			int v; iss >> v;
			if(v == 0) break;
			v--;
			int u;
			while(iss >> u) {
				u--;
				addEdge(u, v);
			}
		}
		findAps(n);
		int cont = 0;
		for(int i = 0; i < n; i++) cont += ap[i];
		cout << cont << endl;
		clearGraph(n);
	}
}