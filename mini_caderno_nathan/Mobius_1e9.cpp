#pragma GCC optimize("Ofast")
#pragma GCC optimize ("unroll-loops")
// #pragma GCC target("sse,sse2,sse3,ssse3,sse4,popcnt,abm,mmx,avx,tune=native")

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

const int N = 5e6 + 10;

struct prefix_mul {

    typedef long long (*func) (long long);

    func p_f, p_g, p_c;
    long long n, th, inv;
    std::unordered_map <long long, long long> mem;

    prefix_mul (func p_f, func p_g, func p_c) : p_f (p_f), p_g (p_g), p_c (p_c) {}

    long long calc (long long x) {
        if (x <= th) return p_f (x);
        auto d = mem.find (x);
        if (d != mem.end ()) return d -> second;
        long long ans = 0;
        for (long long i = 2, la; i <= x; i = la + 1) {
            la = x / (x / i);
            ans = ans + (p_g (la) - p_g (i - 1)) * calc (x / i);
        }
        ans = p_c (x) - ans; ans = ans / inv;
        return mem[x] = ans;
    }

    long long solve (long long n, long long th) {
        if (n <= 0) return 0;
        prefix_mul::n = n; prefix_mul::th = th;
        inv = p_g (1);
        return calc (n); 
    }

};

const int mod = 998244353;

int fast_exp(int b, int e) {
   if(e == 0) return 1;
   
   int ans = fast_exp(b, e / 2);
   ans = (ans * ans) % mod;
   if(e % 2 == 1) ans = (ans * b) % mod;
   
   return ans;
}

int n, U[N];

void mobius(){
    for(int i = 2; i < N; i++) U[i] = 2;
    U[1] = 1;
    for(int i = 2; i < N; i++){
        if(U[i] == 2)
            for(int j = i; j < N; j += i)
                if(U[j]){
                    if(U[j] == 2) U[j] = 1;
                    U[j] *= j / i % i ? -1 : 0;
                }
    }
    for(int i = 1; i < N - 1; i++) U[i + 1] += U[i];
}

int P_g(int x){
    return x;
}

int P_c(int x){
    return 1;
}

int P_f(int x){
    assert(x < N);
    return U[x];
}

unordered_map<int, int> memo;

int pref(int x, prefix_mul &Sf){
    if(memo[x]) return memo[x];
    int ans = 0;
    for (int i = 1, la; i <= x; i = la + 1) {
        la = x / (x / i);
        ans = ans + (Sf.solve(la, N) - Sf.solve(i - 1, N)) * (x / i) * (x / i);
    }
    return ans;
}

int sum_phi(int l, int r, prefix_mul &Sf){
    int dir = pref(r, Sf);

    l--;
    int esq = pref(l, Sf);
    dir -= esq;
    return dir;
}

void solve(){
    prefix_mul Sf(P_f, P_g, P_c);

    cin >> n;

    int quant = sum_phi(1, n, Sf), ways = n;

    for (int i = 2, la; i < n; i = la + 1) {
        la = n / (n / i);
        if(la == n) break;
        int range = sum_phi(i, la, Sf);
        ways = ways * fast_exp(n / i, range) % mod;
    }

    cout << quant % mod << " " << ways << "\n";
}

signed main(){
    IOS;
    mobius();
    solve();
}