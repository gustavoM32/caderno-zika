/*
Latam 2014 - J
https://www.beecrowd.com.br/judge/en/problems/view/1752
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

mt19937 rng((int) chrono::steady_clock::now().time_since_epoch().count());

const int mod = 1e9 + 7;
const int MAXN = 1010;

struct Node{
    int peso, x, y;

    Node(){}

    Node(int PESO, int X, int Y){
        peso = PESO;
        x = X;
        y = Y;
    }

    Node operator +(const Node &o) const{
        Node ans;
        if(peso < o.peso) ans = Node(peso, x, y);
        else ans = Node(o.peso, o.x, o.y);
        return ans;
    }
};

int n,m,Q;
Node a[MAXN][MAXN],st[2*MAXN][2*MAXN];

Node op(Node a, Node b){
    return a + b;
}

void build(){
    for(int i = 0; i < n; i++) for(int j = 0; j < m; j++)st[i+n][j+m]=a[i][j];
    for(int i = 0; i < n; i++) for(int j=m-1;j;--j)
        st[i+n][j]=op(st[i+n][j<<1],st[i+n][j<<1|1]);
    for(int i=n-1;i;--i) for(int j = 0; j < 2 * m; j++)
        st[i][j]=op(st[i<<1][j],st[i<<1|1][j]);
}
void upd(int x, int y, Node v){
    st[x+n][y+m]=v;
    for(int j=y+m;j>1;j>>=1)st[x+n][j>>1]=op(st[x+n][j],st[x+n][j^1]);
    for(int i=x+n;i>1;i>>=1)for(int j=y+m;j;j>>=1)
        st[i>>1][j]=op(st[i][j],st[i^1][j]);
}

//essa query vai de x0, y0 ate x1 - 1, y1 - 1 !!!
Node query(int x0, int x1, int y0, int y1){
    Node r = Node(infinite, 0, 0); // definir elemento neutro da query!!!
    for(int i0=x0+n,i1=x1+n;i0<i1;i0>>=1,i1>>=1){
        int t[4],q=0;
        if(i0&1)t[q++]=i0++;
        if(i1&1)t[q++]=--i1;
        for(int k = 0; k < q; k++) for(int j0=y0+m,j1=y1+m;j0<j1;j0>>=1,j1>>=1){
            if(j0&1)r=op(r,st[t[k]][j0++]);
            if(j1&1)r=op(r,st[t[k]][--j1]);
        }
    }
    return r;
}

int casa[MAXN][MAXN][3];
// 0 -> cost
// 1 -> vertical
// 2 -> horizontal

int operacao(pair<int, int> ini, pair<int, int> fim){
    build();
    vector dist(n, vector<int>(m, infinite));

    dist[ini.first][ini.second] = 0;
    upd(ini.first, ini.second, Node(infinite, 0, 0));

    set< array<int, 3> > s;
    s.insert({0, ini.first, ini.second});

    while(s.size()){
        auto it = *s.begin(); s.erase(it);

        int x = it[1], y = it[2];
        if(it[0] > dist[x][y]) continue;

        int R = casa[x][y][1], C = casa[x][y][2];
        while(true){
            Node it1 = query(max(x - R, 0ll), min(x + R + 1, n), max(y - C, 0ll), min(y + C + 1, m));
            if(it1.peso > infinite / 2) break;

            dist[it1.x][it1.y] = dist[x][y] + it1.peso;
            s.insert({dist[it1.x][it1.y], it1.x, it1.y});
            upd(it1.x, it1.y, Node(infinite, 0, 0));
        }
    }
    int ans = dist[fim.first][fim.second] + casa[ini.first][ini.second][0] - casa[fim.first][fim.second][0];
    if(ans > infinite / 2) return -1;
    return ans;
}

void solve(){
    cin >> n >> m >> Q;

    for(int k = 0; k < 3; k++)
        for(int i = 0; i < n; i++)
            for(int j = 0; j < m; j++)
                cin >> casa[i][j][k];

    for(int i = 0; i < n; i++)
        for(int j = 0; j < m; j++)
            a[i][j] = Node(casa[i][j][0], i, j);

    vector< pair<int, int> > province;
    for(int i = 0; i < Q; i++){
        int x, y; cin >> x >> y; 
        x--; y--;
        province.pb(mp(x, y));
    }

    for(int i = 0; i < Q - 2; i++) cout << operacao(province[i], province[i + 1]) << " ";
    cout << operacao(province[Q - 2], province[Q - 1]) << "\n";
}

signed main(){
    IOS;
    solve();
}