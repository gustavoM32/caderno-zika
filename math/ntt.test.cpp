#define PROBLEM "https://codeforces.com/contest/1613/problem/F"
#define IGNORE // incompatible judge
/**********/
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
/**********/
const int MOD = 998244353;

namespace NTT {
    const long long mod = 998244353;
    const long long root = 15311432;
    const long long root_1 = 469870224;
    const long long root_pw = 1 << 23;
 
    long long fastxp(long long n, long long e){
        long long ans = 1, pwr = n;
        while(e){
            if(e%2)  ans = ans * pwr % MOD;
            e /= 2;
            pwr = pwr * pwr % MOD;
        }
        return ans % MOD;
    }
 
 
    void fft(vector<long long> & a, bool invert) {
        long long n = a.size();
 
        for (long long i = 1, j = 0; i < n; i++) {
            long long bit = n >> 1;
            for (; j & bit; bit >>= 1)
                j ^= bit;
            j ^= bit;
 
            if (i < j) swap(a[i], a[j]);
        }
 
        for (long long len = 2; len <= n; len <<= 1) {
            long long wlen = invert ? root_1 : root;
            for (long long i = len; i < root_pw; i <<= 1)
                wlen = (long long)(1LL * wlen * wlen % mod);
 
            for (long long i = 0; i < n; i += len) {
                long long w = 1;
                for (long long j = 0; j < len / 2; j++) {
                    long long u = a[i+j], v = a[i+j+len/2] * w % mod;
                    a[i+j] = u + v < mod ? u + v : u + v - mod;
                    a[i+j+len/2] = u - v >= 0 ? u - v : u - v + mod;
                    w = (long long)(1LL * w * wlen % mod);
                }
            }
        }
 
        if (invert) {
            long long n_1 = fastxp(n, mod - 2);
            for (long long & x : a) x = x * n_1 % mod;
        }
    }
 
    vector<long long> multiply(vector<long long> &a, vector<long long> &b) {
        vector<long long> fa(a.begin(), a.end()), fb(b.begin(), b.end());
        long long sz = a.size() + b.size() - 1, n = 1;
        while (n < sz) n <<= 1;
 
        fa.resize(n), fb.resize(n);
        fft(fa, 0), fft(fb, 0);
        for (long long i = 0; i < fa.size(); i++) fa[i] = fa[i] * fb[i] % mod;
 
        fft(fa, 1);
        fa.resize(sz);
        return fa;
    }
};

int n, degree[300010];

vector<int> divide_conquer(int ini, int fim){
    if(ini == fim) return {1, degree[ini]};
    
    int meio = (ini + fim) / 2;
    auto esq = divide_conquer(ini, meio);
    auto dir = divide_conquer(meio + 1, fim);

    return NTT::multiply(esq, dir);
}

void solve(){
    cin >> n;
    for(int i = 1; i < n; i++){
        int a, b; cin >> a >> b;
        degree[a]++;
        degree[b]++;
    }
    for(int i = 2; i <= n; i++) degree[i]--;

    // queremos fazer Produtorio (1 + x * degree[i])
    vector<int> fat(n + 10);
    fat[0] = fat[1] = 1;
    for(int i = 2; i <= n; i++) fat[i] = fat[i - 1] * i % MOD;

    auto a = divide_conquer(1, n);
    int ans = 0;
    for(int i = 0, sinal = 1; i <= n; i++, sinal *= -1){
        ans = (ans + sinal * a[i] * fat[n - i] % MOD + MOD) % MOD;
    }
    cout << ans << "\n";
}

signed main(){
    IOS;
    solve();
}