#define IGNORE // incompatible judge
/**********/
#include <bits/stdc++.h>
using namespace std;

#define FASTIO ios_base::sync_with_stdio(false);cin.tie(NULL)
#define pb push_back
#define mp make_pair
#define sz(x) int(x.size())
#define trace(x) cerr << #x << ": " << x <<endl;
/**********/
typedef long long ll;
const ll N = 512;
const ll INF = 1ll << 61;
const ll MOD = 1e9 + 7;

#include "gauss_mod2.cpp"

int main() {
  int n, m;
  cin >> n >> m;
  vector<bitset<N>> a(n, bitset<N>(0));

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; j++) {
      int x;
      cin >> x;
      if (x) a[i].set(j);
    }
  }

  for (int i = 0; i < n; i++) {
    int x;
    cin >> x;
    if (x) a[i].set(m);
  }

  Gauss<N> g(n, m, a);
  bitset<N> ans;
  int s = g.solve(ans);

  if (s == 0) cout << "no solution\n";
  else {
    if (s == 1) cout << "unique solution\n";
    else cout << "multiple solutions\n";

    for (int i = 0; i < m; i++) {
      cout << ans[i] << " ";
    }
    cout << "\n";
  }

  return 0;
}
