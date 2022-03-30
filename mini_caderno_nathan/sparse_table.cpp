#include <bits/stdc++.h>
using namespace std;

#define pb push_back
#define mp make_pair

#define int long long
#define IOS ios_base::sync_with_stdio(0);cin.tie(0)
#define PRECISION cout.precision(3); cout.setf(ios::fixed);
#define fr(i,n) for(int i = 0; i<n; i++)
#define sz(v) (int)(v.size())

#define infinite 112345678912345
#define db cout << "Debug" << "\n";
#define dbg(x)  cout << #x << " = " << x << "\n"

int log_floor(int n){
    if(n == 0) return 1;
    return 31-__builtin_clz(n);
}

int oper(int a, int b){
    return __gcd(a, b);
}

struct sparse_table{
    int exp2;
    int n;
    vector<vector<int>> mat;
    sparse_table(){}
    sparse_table(vector<int> v){
        n = sz(v);
        exp2 = log_floor(n)+1;
        mat.resize(exp2);
        mat[0].resize(n);
        fr(i,n) mat[0][i] = v[i];
        for(int k = 1; k<exp2; k++){
            mat[k].resize(n);
            for(int i = 0; i+(1<<k)<=n; i++){
                mat[k][i] = oper(mat[k-1][i],mat[k-1][i+(1<<(k-1))]);
            }
        }
    }
    //query fechada [l,r]
    int qry(int l, int r){
        assert(l<=r and l>=0 and r<n);
        int k = log_floor(r-l+1);
        return oper(mat[k][l],mat[k][r-(1<<k)+1]);
    }
}; //end sparse_table

int n;

void solve(){
    cin >> n;
    vector<int> a;
    for(int i = 0; i < n; i++){
        int x; cin >> x;
        a.pb(x);
    }

    sparse_table ST(a);

    vector<int> ans(n + 10);
    int esq = 0, cur = 0;
    for(int i = 0; i < n; i++){
        bool deu_ruim = false;
        
        int d = a[i], dir = i;
        while(dir >= esq){
            int ini = esq, fim = dir;
            while(ini <= fim){
                int meio = (ini + fim) / 2;
                if(ST.qry(meio, i) < d) ini = meio + 1;
                else fim = meio - 1;
            }

            int l = i - dir + 1, r = i - ini + 1;

            if(l <= d and d <= r){
                deu_ruim = true;
                break;
            }
            dir = ini - 1;
            if(ini) d = ST.qry(ini - 1, i);
        }

        if(i) ans[i] = ans[i - 1];
        if(deu_ruim){
            ans[i]++;
            esq = i + 1;
        }
    }

    for(int i = 0; i < n; i++) cout << ans[i] << " ";
}

signed main(){
    IOS;
    solve();
}