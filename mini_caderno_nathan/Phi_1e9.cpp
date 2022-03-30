/*
https://matcomgrader.com/problem/9592/straight-lines/
*/
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

const int N = 1e6;
const int MOD = 998244353;
const int MAXN = 5e6;

bool is_prime[MAXN];
int cnt,phi[MAXN],prime[MAXN];
int f[MAXN];
map<int,int> mp;


int mul_mod(int a,int i)
{
    int s=0;a%=MOD;
    while(i)
    {
        if(i&1) s=(s+a)%MOD;
        a=(a+a)%MOD;
        i>>=1;
    }
    return s;
}
 
int pow_mod(int a,int i)
{
    int s=1;
    while(i)
    {
        if(i&1) s=mul_mod(s,a);
        a=mul_mod(a,a);
        i>>=1;
    }
    return s;
}
void genphi(int n)
{
    int p=0;
    memset(phi,0,sizeof(phi));
    phi[1]=1;
     for(int i=0;i<=n;i++) is_prime[i]=true;
    is_prime[0]=is_prime[1]=false;
    for(int i=2;i<=n;i++)
    {
        if(is_prime[i]) {prime[p++]=i; phi[i]=i-1;}
        for(int j=0;j<p;j++)
        {
            if(prime[j]*i>n) break;
            is_prime[prime[j]*i]=false;
            phi[i*prime[j]]=phi[i]*(i%prime[j]?prime[j]-1:prime[j]);
            if(i%prime[j]==0) break;
        }
    }
    // for(int i = 1; i <= 10; i++) cout << phi[i] << " ";
    for(int i=1;i<=n;i++) f[i]=(f[i-1]+phi[i]);
}
int calc(int x)
{
    if(x<=N) return f[x];
    if(mp.find(x)!=mp.end()) return mp[x];
    int ans = x * (x + 1) / 2;
    for(int i=2,r;i<=x;i=r+1)
    {
        r=x/(x/i);
        ans -= calc(x/i)*(r-i+1);
    }
    return mp[x]=ans;
}

int sum_phi(int l, int r){
    return 2 * (calc(r) - calc(l - 1)) - (l == 1 ? 1 : 0);
}

pair<int, int> solve(int n){
    int quant = sum_phi(1, n), ways = n;

    for (int i = 2, la; i < n; i = la + 1) {
        la = n / (n / i);
        if(la == n) break;
        int range = sum_phi(i, la);
        ways = mul_mod(ways, pow_mod(n / i, range)) % MOD;
    }

    return make_pair(quant, ways);
}

pair<int, int> brute(int n){
    int quant = 0, ways = 1;
    for(int i = 1; i <= n; i++){
        for(int j = 1; j <= n; j++){
            if(__gcd(i, j) == 1){
                quant++;
                ways = ways * (n / max(i, j)) % MOD;
            }
        }
    }
    return make_pair(quant, ways);
}

int n;

signed main()
{
    genphi(N);
    cin >> n;

    auto it = solve(n);
    cout << it.first << " " << it.second << "\n";
}