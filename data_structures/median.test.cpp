#include <bits/stdc++.h>
#define PROBLEM "https://www.hackerrank.com/challenges/find-the-running-median/problem"
using namespace std;

#define fastio ios_base::sync_with_stdio(0);cin.tie(0)
#define pb push_back
#define mp make_pair
#define sz(x) int(x.size())
#define trace(x) cerr << #x << ": " << x <<endl;

typedef long long ll;

const ll N = 1e6;
const ll INF = 1LL << 61;
const ll MOD = 1e9 + 7;

#include "median.cpp"

int main() {
    Median m;
    int n;
    cin >> n;
    vector<int> ans;
    for (int i = 0; i < n; i++) {
        ll a;
        cin >> a;
        m.insert(a);
        cout << setprecision(1) << fixed << m.median() << "\n";
    }

    return 0;
}