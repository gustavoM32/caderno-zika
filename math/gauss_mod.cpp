/* Solves a system of linear equations in Z_mod
 *
 * Complexity: O(n^3)
 */
ll gcdExt(ll a, ll b, ll &x, ll &y) {
	if (b == 0) {
		x = 1;
		y = 0;
		return a;
	}
	ll x1, y1;
	ll res = gcdExt(b, a % b, x1, y1);
	x = y1;
	y = x1 - y1 * (a / b);
	return res;
}

ll inv(ll a, ll m = MOD) {
	ll x, y;
	ll g = gcdExt(a, m, x, y);
	if (g != 1) return -1;
	return (m + x % m) % m;
}

struct Gauss {
  int n, m;
  vector<int> pos;
  int rank = 0;
  ll mod;
  vector<vector<ll>> a;

  // n equations, m-1 variables, last column is for coefficients
  Gauss(int n, int m, ll mod, vector<vector<ll>> &a) : n(n), m(m), mod(mod), a(a) {
    pos.assign(m, -1);
  }

  /* if a solution exists, it will be stored in ans 
    0 - no solution
    1 - unique solution
    2 - infinite number of solutions */
  int solve(vector<ll> &ans) {
    for (int col = 0, row = 0; col < m && row < n; col++) {
      int sel = row;
      for (int i = row+1; i < n; i++) {
        if (a[i][col] > 0) {
          sel = i;
          break;
        }
      }

      if (a[sel][col] == 0) continue;

      swap(a[sel], a[row]);

      pos[col] = row;

      for (int i = 0; i < n; i++) {
        if (i != row) {
          ll mult = a[i][col] * inv(a[row][col], mod) % mod;
          a[i][col] = 0;
          for (int j = col+1; j <= m; j++) {
            a[i][j] = (a[i][j] + mod - a[row][j] * mult % mod) % mod;
          }
        }
      }

      ++row, ++rank;
    }

    ans.assign(m, 0);

    bool multiple = false;

    for (int i = 0; i < m; i++) {
      if (pos[i] != -1) ans[i] = a[pos[i]][m] * inv(a[pos[i]][i], mod) % mod;
      else multiple = true;
    }

    for (int i = 0; i < n; i++) {
      ll sum = 0.0;
      for (int j = 0; j < m; j++) {
        sum = (sum + ans[j] * a[i][j] % mod) % mod;
      }
      if (sum != a[i][m]) return 0;
    }

    if (multiple) return 2;

    return 1;
  }
};
