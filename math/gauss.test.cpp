#define IGNORE
#include <bits/stdc++.h>
using namespace std;

#define FASTIO ios_base::sync_with_stdio(false);cin.tie(NULL)
#define pb push_back
#define mp make_pair
#define sz(x) int(x.size())
#define trace(x) cerr << #x << ": " << x <<endl;

typedef long long ll;
const ll N = 512;
const ll INF = 1ll << 61;
const ll MOD = 1e9 + 7;

#include "gauss.cpp"

int main() {
  int n, m;
  cin >> n >> m;
  vector<vector<double>> a(n, vector<double>(m+1));
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; j++) {
      cin >> a[i][j];
    }
  }

  for (int i = 0; i < n; i++) cin >> a[i][m];
    
  Gauss g(n, m, a);

  vector<double> ans(n);
  g.solve(ans);

  for (int i = 0; i < n; i++) {
    cout << ans[i] << " ";
  }
  cout << "\n";

  return 0;
}
