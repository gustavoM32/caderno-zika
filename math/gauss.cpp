/* Solves a system of linear equations
 *
 * Complexity: O(n^3)
 */
struct Gauss {
  int n, m;
  vector<int> pos;
  int rank = 0;
  vector<vector<double>> a;
  const double EPS = 1e-9;

  Gauss(int n, int m, vector<vector<double>> &a) : n(n), m(m), a(a) {
    pos.assign(m, -1);
  }

  /* if a solution exists, it will be stored in ans 
    0 - no solution
    1 - unique solution
    2 - infinite number of solutions */
  int solve(vector<double> &ans) {
    for (int col = 0, row = 0; col < m && row < n; col++) {
      int sel = row;
      for (int i = row+1; i < n; i++) {
        if (abs(a[i][col]) > abs(a[sel][col])) sel = i;
      }

      if (abs(a[sel][col]) < EPS) continue;

      swap(a[sel], a[row]);

      pos[col] = row;

      for (int i = 0; i < n; i++)
        if (i != row) {
          double mult = a[i][col] / a[row][col];
          a[i][col] = 0.0;
          for (int j = col+1; j <= m; j++) {
            a[i][j] -= a[row][j] * mult;
          }
        }

      ++row, ++rank;
    }

    ans.assign(m, 0);

    bool multiple = false;

    for (int i = 0; i < m; i++) {
      if (pos[i] != -1) ans[i] = a[pos[i]][m] / a[pos[i]][i];
      else multiple = true;
    }

    for (int i = 0; i < n; i++) {
      double sum = 0.0;
      for (int j = 0; j < m; j++) {
        sum += ans[j] * a[i][j];
      }
      if (abs(sum - a[i][m]) > EPS) return 0;
    }

    if (multiple) return 2;

    return 1;
  }
};
