#define PROBLEM "https://onlinejudge.org/index.php?option=com_onlinejudge&Itemid=8&category=16&page=show_problem&problem=1369"
#define IGNORE // incompatible judge
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
const double EPS = 1e-7;
const ll INF = 1LL << 61;
const ll MOD = 1e9 + 7;

#include "roots_newton.cpp"

double evalPoly(vector<double> &p, double x) {
    double y = 0.0;
    double xi = 1;
    for (double pi : p) {
        y += pi * xi;
        xi *= x;
    }
    return y;
}

void divideByRoot(vector<double> p, vector<double> &q, double x) {
    int n = p.size();
    q.resize(n-1);
    for (int i = n - 1; i > 0; i--) {
        q[i-1] = p[i];
        p[i-1] += q[i-1] * x;
    }
}

void derivate(vector<double> &p, vector<double> &dp) {
    int n = p.size();
    dp.resize(n-1);
    for (int i = 1; i < n; i++) {
        dp[i-1] = i * p[i];
    }
}

vector<double> p;
vector<double> dp;

double f(double x) {
    return evalPoly(p, x);
}

double df(double x) {
    return evalPoly(dp, x);
}

int main() {
    int n;
    cin >> n;
    int cc = 1;

    while (n > 0) {
        p.resize(n+1);
        for (int i = n; i >= 0; i--) {
            cin >> p[i];
        }
        vector<double> roots;
        for (int i = 0; i < n; i++) {
            derivate(p, dp);
            double r = findRoot();
            roots.push_back(r);
            divideByRoot(p, p, r);
        }
        sort(roots.begin(), roots.end());

        cout << "Equation " << cc++ << ":";
        for (int i = 0; i < roots.size(); i++) {
            cout << " " << setprecision(4) << fixed << roots[i];
        }
        cout << "\n";
        cin >> n;
    }

    return 0;
}

