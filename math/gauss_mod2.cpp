/* Solves a system of linear equations in Z_2, it is faster than gauss_mod
 * by the use of bitset.
 *
 * Complexity: O(n^3)
 */
template<int M>
struct Gauss {
  int n, m;
  array<int, M> pos;
  int rank = 0;
  vector<bitset<M>> a;

  Gauss(int n, int m, vector<bitset<M>> &a) : n(n), m(m), a(a) {
    pos.fill(-1);
  }

  int solve(bitset<N> &ans) {
    for (int col = 0, row = 0; col < m && row < n; col++) {
      int one = -1;
      for (int i = row; i < n; i++) {
        if (a[i][col]) {
          one = i;
          break;
        }
      }

      if (one == -1) { continue; }

      swap(a[one], a[row]);

      pos[col] = row;

      for (int i = row + 1; i < n; i++)
        if (a[i][col])
          a[i] ^= a[row];

      ++row, ++rank;
    }

    ans.reset();

    for (int i = m - 1; i >= 0; i--) {
      if (pos[i] == -1) ans[i] = true;
      else {
        int k = pos[i];
        for (int j = i + 1; j < m; j++) if (a[k][j]) ans[i] = ans[i] ^ ans[j];
        ans[i] = ans[i] ^ a[k][m];
      }
    }

    for (int i = rank; i < n; i++) if (a[i][m]) return 0; 

    return 1;
  }
};
