/* Solves any flow problem with demands (also known as lower bounds) in the flow of
 * some edges.
 * It only does a graph transformation and applies the flow algorithm on it, thus this
 * template depends of another flow algorithm template. In this example, that algorithm
 * is minimum-cost maximum-flow, but that's not always necessary. Adapting this template
 * to other flow algorithms should be easy.
 *
 * Constructor:
 * demands(n, s, t)
 * n - number of nodes in the original flow network.
 * s - source of the original flow network.
 * t - sink of the original flow network.
 *
 * Methods:
 * - addEdge(u, v, cap, dem, cost)
 *   adds an edge to the original flow network from u to v with capacity `cap`, flow
 *   demand `dem` and cost `cost`. For problems without cost, the last parameter should
 *   be erased.
 *
 * - getFlow()
 *   finishes building the auxiliary graph and returns the result of applying the flow
 *   algorithm on it.
 *
 * Complexity: depends on the flow algorithm used, but take into account that this method
 *             adds O(V) edges to the network flow.
 */

struct demands {
	vector<ll> din, dout;
	mcf g;
	int n, s, t;
	ll base = 0;
	demands(int n, int _s, int _t) : n(n), s(n), t(n + 1), din(n, 0), dout(n, 0) {
		g = mcf(n + 2, s, t);
		g.addEdge(_t, _s, INF, 0);
	};
	void addEdge(int u, int v, ll cap, ll dem, ll cost) {
		din[v] += dem;
		dout[u] += dem;
		base += dem * cost;
		g.addEdge(u, v, cap - dem, cost);
	}
	pair<ll, ll> getFlow() {
		for(int i = 0; i < n; i++) {
			if(din[i] > 0) g.addEdge(s, i, din[i], 0);
			if(dout[i] > 0) g.addEdge(i, t, dout[i], 0);
		}
		auto ans = g.getFlow();
		return {ans.first, base + ans.second};
	}
};