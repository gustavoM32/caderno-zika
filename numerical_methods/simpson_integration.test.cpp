#define PROBLEM "https://www.beecrowd.com.br/judge/en/problems/view/1297"
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

const ll N = 1e3;
const double EPS = 1e-7;
const ll INF = 1LL << 61;
const ll MOD = 1e9 + 7;

#include "simpson_integration.cpp"

vector<int> p1, p2, q1, q2;
double d;
int W, D, A, K;

double poly(vector<int> p, double x) {
    double y = 0;
    double xi = 1;
    for (int pi : p) {
        y += xi * pi;
        xi *= x;
    }
    return y;
}

double f(double x) {
    double y1 = poly(p1, x) / poly(q1, x);
    double y2 = poly(p2, x) / poly(q2, x);
    y1 = max(y1, -d);
    y2 = max(y2, -d);
    return y1 - y2;
}

double bb(double l, double r) {
    while (r - l > EPS) {
        d = (l + r) / 2;
        double dd = integral(0, W, N);
        if (integral(0, W, N) < A) {
            l = d;
        } else {
            r = d;
        }
    }
    return l;
}

void read(vector<int> &p) {
    p.resize(K+1);
    for (int i = 0; i <= K; i++) {
        cin >> p[i];
    }
}

int main() {
    while (cin >> W >> D >> A >> K) {
        read(p1);
        read(q1);
        read(p2);
        read(q2);
        cout << setprecision(5) << fixed << bb(0, D) << "\n";
    }
}

