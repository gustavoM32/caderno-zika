/* Hungarian method algorithm that solves the bipartite perfect matching of minimum cost.
 * Can solve the problem of maximum cost with minor changes.
 *
 * Constructor:
 * hungarian(n, m)
 * n - (input) size of color class X
 * m - (input) size of color class Y
 *
 * Methods:
 * - set(x, y, c)
 *   sets the cost c for the edge between x \in X to y \in Y
 * - assign()
 *   returns the cost of an optimal matching and fills the vectors matchl and matchr
 *   with the assignment of such matching 
 * Complexity: O(V^3) where V is the number of vertex of the graph
 */

typedef long double ld;
const ld INF = 1e100; // for maximization set INF to 0 and negate costs

bool zero(ld x) {
	return fabs(x) < 1e-9; // change to x == 0 for integer types
}

struct hungarian {
	int n;
	vector<vector<ld> > cs;
	vector<int> matchl, matchr;
	hungarian(int _n, int _m) : n(max(_n, _m)), cs(n, vector<ld>(n)), matchl(n), matchr(n) {
		for(int x = 0; x < _n; x++)
			for(int y = 0; y < _m; y++)
				cs[x][y] = INF;
	}
	void set(int x, int y, ld c) {
		cs[x][y] = c;
	}
	ld assign() {
		int mat = 0;
		vector<ld> ds(n), y(n), z(n);
		vector<int> dad(n), vis(n);
		for(int i = 0; i < n; i++) {
			matchl[i] = matchr[i] = -1;
			y[i] = *min_element(cs[i].begin(), cs[i].end());
		}
		for(int j = 0; j < n; j++) {
			z[j] = cs[0][j] - y[0];
			for(int i = 1; i < n; i++)
				z[j] = min(z[j], cs[i][j] - y[i]);
		}
		for(int i = 0; i < n; i++)
			for(int j = 0; j < n; j++)
				if(matchr[j] == -1 && zero(cs[i][j] - y[i] - z[j])) {
					matchl[i] = j;
					matchr[j] = i;
					mat++;
					break;
				}
		for(;mat < n; mat++) {
			int s = 0, j, i;
			while(matchl[s] != -1) s++;
			for(int i = 0; i < n; i++) {
				dad[i] = -1;
				vis[i] = 0;
			}
			for(int k = 0; k < n; k++)
				ds[k] = cs[s][k] - y[s] - z[k];
			while(1) {
				j = -1;
				for(int k = 0; k < n; k++)
					if(!vis[k] && (j == -1 || ds[k] < ds[j]))
						j = k;
				vis[j] = 1;
				i = matchr[j];
				if(i == -1)
					break;
				for(int k = 0; k < n; k++) if(!vis[k]) {
					auto new_ds = ds[j] + cs[i][k] - y[i] - z[k];
					if(ds[k] > new_ds) {
						ds[k] = new_ds;
						dad[k] = j;
					}
				}
			}
			for(int k = 0; k < n; k++) if(k != j && vis[k]) {
				auto w = ds[k] - ds[j];
				z[k] += w;
				y[matchr[k]] -= w;
			}
			y[s] += ds[j];
			while(dad[j] != -1) {
				matchl[matchr[j] = matchr[dad[j]]] = j;
				j = dad[j];
			}
			matchr[j] = s;
			matchl[s] = j;
		}
		ld value = 0;
		for(int i = 0; i < n; i++)
			value += cs[i][matchl[i]];
		return value;
	}
};