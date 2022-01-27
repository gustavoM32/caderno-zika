// https://codeforces.com/contest/632/problem/E
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


namespace FFT {
    using cd = complex<double>;
 
    #define PI acos(-1)
    void fft(vector<cd> &a, int logN, int sign) {
        for (int i = 0; i < a.size(); i++) {
            int mask = 0;
            for (int j = 0; j < logN; j++) {
                if (i >> j & 1) mask |= (1 << logN - j - 1);
            }
            if (i < mask) swap(a[i], a[mask]);
        }
 
        for (int k = 1; k < a.size(); k <<= 1) {
            cd wDelta = polar(1.0, PI * sign / k);
            cd w(1, 0);
            for (int i = 0; i < k; i++, w *= wDelta) {
                for (int j = i; j < a.size(); j += k * 2) {
                    cd foo = a[j];
                    cd bar = w * a[j + k];
                    a[j] = foo + bar;
                    a[j + k] = foo - bar;
                }
            }
        }
    }
 
    vector<long long> multiply(vector<long long> &a, vector<long long> &b) {
        int sz = a.size() + b.size() - 1;
        int logN = 0;
        while ((1 << logN) < sz) logN++;
 
        vector<cd> aa(1 << logN, 0), bb(1 << logN, 0);
 
        for (int i = 0; i < a.size(); i++) aa[i] = a[i];
        for (int i = 0; i < b.size(); i++) bb[i] = b[i];
 
        fft(aa, logN, 1);
        fft(bb, logN, 1);
 
        for (int i = 0; i < aa.size(); i++) {
            aa[i] *= bb[i];
        }
 
        fft(aa, logN, -1);
        vector<long long> res(sz);
        for (int i = 0; i < sz; i++) {
            res[i] = (long long)(aa[i].real() / (1 << logN) + 0.5);
        }
        return res;
    }
    
    vector<long long> fastxp(vector<long long> a, long long exp){
        if(exp == 1) return a;
 
        vector<long long> aux = fastxp(a, exp/2);
 
        aux = multiply(aux, aux);
 
        if(aux.size() > 1000010) 
            aux.resize(1000010);
 
        if(exp%2 == 1)
            aux = multiply(aux, a);
 
        if(aux.size() > 1000010) 
            aux.resize(1000010);
        
        // evitar overflow ou erro de precisao
        for(int i = 0; i < aux.size(); i++) if(aux[i]) aux[i] = 1;
 
        return aux;
    }
};

int n, k, a[1010];

void solve(){
    cin >> n >> k;
    vector<int> rec(1001);
    for(int i = 1; i <= n; i++){
        cin >> a[i];
        rec[a[i]] = 1;
    }

    vector<int> ans = FFT::fastxp(rec, k);
    for(int i = 0; i < ans.size(); i++) if(ans[i] != 0) cout << i << " ";
}

signed main(){
    IOS;
    solve();
}