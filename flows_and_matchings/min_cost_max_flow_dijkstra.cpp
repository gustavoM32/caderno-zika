/* Solves the minimum-cost maximum-flow problem.
 *
 * Constructor:
 * mcf(n, s, t)
 * n - number of nodes in the flow graph.
 * s - source of the flow graph.
 * t - sink of the flow graph.
 *
 * Methods:
 * - addEdge(u, v, cap, cost)
 *   adds a directed edge from u to v with capacity `cap` and cost `cost`.
 * - getFlow()
 *   returns a pair of integers in which the first value is the maximum flow and the
 *   second is the minimum cost to achieve this flow.
 *
 * Complexity: There are two upper bounds to the time complexity of getFlow
 *              - O(max_flow * (E log V))
 *              - O(V * E * (E log V))
 */

const ll INF = 1LL << 60;

struct mcf {
	int n, s, t;
	ll cost, fl;
	vector<int> first, prev;
	vector<ll> dist;
	struct edge {
		int to, next;
		ll cap, cost;
		edge(int _to, ll _cap, ll _cost, int _next)
            : to(_to), cap(_cap), cost(_cost), next(_next) {};
	};
	vector<edge> g;
	mcf() {}
	mcf(int _n, int _s,int _t) : n(_n), s(_s), t(_t), fl(0), cost(0) {
		first.resize(n, -1);
		dist.resize(n);
		prev.resize(n);
		g.reserve(n*n);
	};
	void addEdge(int u, int v, ll cap, ll cost) {
		g.pb(edge(v, cap, cost, first[u]));
		first[u] = sz(g) - 1;
		g.pb(edge(u, 0, -cost, first[v]));
		first[v] = sz(g) - 1;
	}
	bool augment() {
		dist.assign(n, INF);
		dist[s] = 0;
		priority_queue<pair<ll, int> > q;
		q.push({0, s});
		while(!q.empty()) {
			if(dist[t] < INF) break;
			ll d, u;
			tie(d, u) = q.top(); 
			d *= -1;
			q.pop();
			if(dist[u] < d) continue;
			for(int e = first[u]; e != -1; e = g[e].next) {
				int v = g[e].to;
				ll ndist = d + g[e].cost;
				if(g[e].cap > 0 && dist[v] > ndist) {
					dist[v] = ndist;
					q.push({-ndist, v});
					prev[v] = e;
				}
			}
		}
		return dist[t] < INF;
	}
	pair<ll, ll> getFlow() {
		while(augment()) {
			ll cur = t, curf = INF;
			while(cur != s) {
				int e = prev[cur];
				curf = min(curf, g[e].cap);
				cur = g[e^1].to;
			}
			fl += curf; 
			cost += dist[t] * curf;
			cur = t;
			while(cur != s) {
				int e = prev[cur];
				g[e].cap -= curf;
				g[e^1].cap += curf;
				cur = g[e^1].to;
			}
		}
		return {fl, cost};
	}
};