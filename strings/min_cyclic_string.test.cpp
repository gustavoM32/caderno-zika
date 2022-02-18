#define PROBLEM "https://uva.onlinejudge.org/index.php?option=onlinejudge&page=show_problem&problem=660"
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

#include "min_cyclic_string.cpp"

int main() {
    int n;
    cin >> n;
    while (n--) {
        string s;
        cin >> s;
        cout << min_cyclic_string(s) + 1 << "\n";
    }
    return 0;
}
