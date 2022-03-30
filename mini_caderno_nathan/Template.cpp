#include <bits/stdc++.h>
using namespace std;

#define pb push_back
#define mp make_pair

#define int long long
#define IOS ios_base::sync_with_stdio(0);cin.tie(0)
#define PRECISION cout.precision(3); cout.setf(ios::fixed);

#define infinite 9123456789
#define db cout << "Debug" << "\n";
#define dbg(x)  cout << #x << " = " << x << "\n"

#define ler_entrada freopen("test_input.txt", "r", stdin);

// random prime 9 digits 510857131
const int mod = 1e9 + 7;
const int N = 1e6 + 9;

int fat[N], inv[N];

int fast_exp(int b, int e) {
   if(e == 0) return 1;
   
   int ans = fast_exp(b, e / 2);
   ans = (ans * ans) % mod;
   if(e % 2 == 1) ans = (ans * b) % mod;
   
   return ans;
}

void init() {
   fat[0] = fat[1] = 1;
   inv[0] = inv[1] = 1;
   for(int i = 2; i < N; i++) {
      fat[i] = (i * fat[i - 1])%mod;
      inv[i] = fast_exp(fat[i], mod - 2);
   }
}

int comb(int a, int b) {
   return fat[a] * inv[b] % mod * inv[a - b] % mod;
}

mt19937 rng((int) chrono::steady_clock::now().time_since_epoch().count());
// para pegar random: int x = rng();

// ordered_set
#include <ext/pb_ds/assoc_container.hpp> 
#include <ext/pb_ds/tree_policy.hpp> 
using namespace __gnu_pbds;
 
typedef tree<int, null_type, less<int>, rb_tree_tag, tree_order_statistics_node_update> ordered_set;

ordered_set S;

int getPos(int i) {
    return S.order_of_key(i) + 1;
}
 
void tira(int i) {
    S.erase(i);
}
 
void init(int n) {
    for(int i = 1; i <= n; i++)
        S.insert(i);
}

#pragma GCC optimize("Ofast")
#pragma GCC optimize ("unroll-loops")
#pragma GCC target("sse,sse2,sse3,ssse3,sse4,popcnt,abm,mmx,avx,tune=native")