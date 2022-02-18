#define PROBLEM "https://onlinejudge.org/index.php?option=com_onlinejudge&Itemid=8&page=show_problem&problem=4648"
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

vector<pair<int, int>> edge;
int low[N], disc[N], t_in;
bool bridge[N];
vector<int> adj[N];

void compute(int v, int p=-1) {
	low[v] = disc[v] = ++t_in;
	for(int e : adj[v]) {
		int u = edge[e].first ^ edge[e].second ^ v;
		if(!disc[u]) {
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
	for(int i = 0; i < sz(edge); i++) {
		bridge[i] = 0;
	}
	t_in = 0;
	for(int i = 0; i < n; i++) {
        if(!disc[i]) compute(i);
	}
}

// Must be call between testcases
void clearGraph(int n) {
	for(int i = 0; i < n; i++)
		adj[i].clear();
	edge.clear();
}

int main() {
    fastio;
    int n, m;
    while(cin >> n >> m) {
        if(n == 0) break;
        for(int i = 0; i < m; i++) {
            int u, v; cin >> u >> v;
            addEdge(u, v);
        }
        findBridges(n);
        vector<pair<int, int> > ans;
        for(int i = 0; i < m; i++) if(bridge[i]) {
            ans.pb(edge[i]);
        }
        sort(ans.begin(), ans.end());
        cout << sz(ans);
        for(int i = 0; i < sz(ans); i++) {
			cout << " ";
            cout << ans[i].first << " " << ans[i].second;
		}
		cout << endl;
        clearGraph(n);
    }
}