#include "../base_template.cpp"
#include "angular_sweep.cpp"
#define PROBLEM "https://judge.yosupo.jp/problem/sort_points_by_argument"

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