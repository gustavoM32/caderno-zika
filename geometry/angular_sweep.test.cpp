#define PROBLEM "https://judge.yosupo.jp/problem/sort_points_by_argument"
/**********/
#include <bits/stdc++.h>
using namespace std;

#define fastio ios_base::sync_with_stdio(0);cin.tie(0)
#define pb push_back
#define mp make_pair
#define sz(x) int(x.size())
#define trace(x) cerr << #x << ": " << x <<endl;
/**********/
typedef long long ll;

const ll N = 1e6;
const ll INF = 1LL << 61;
const ll MOD = 1e9 + 7;

#include "angular_sweep.cpp"

point p[N];

int main() {
    fastio;
    int n; cin >> n;
    for(int i = 0; i < n; i++) {
        cin >> p[i].x >> p[i].y;
    }
    sort(p, p + n);
    for(int i = 0; i < n; i++) {
        cout << p[i].x << " " << p[i].y << endl;
    }
}