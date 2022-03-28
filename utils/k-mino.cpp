/* Example code to generate k-minos
 *
 * Complexity: not efficient
 */
int dx[4] = {0, 0, -1, 1};
int dy[4] = {-1, 1, 0, 0};

struct Tab {
    int w, h;
    int tab[10][10];

    Tab() {
        w = 0;
        h = 0;
        for (int i = 0; i < 10; i++) {
            for (int j = 0; j < 10; j++) {
                tab[i][j] = 0;
            }
        }
    }

    Tab(const Tab &t) {
        w = t.w;
        h = t.h;
        for (int i = 0; i < 10; i++) {
            for (int j = 0; j < 10; j++) {
                tab[i][j] = t.tab[i][j];
            }
        }
    }

    void shiftDown() {
        for (int i = h; i > 0; i--) {
            for (int j = 0; j < w; j++) {
                tab[i][j] = tab[i-1][j];
            }
        }
        for (int j = 0; j < w; j++) tab[0][j] = 0;
        h++;
    }

    void shiftRight() {
        for (int j = w; j > 0; j--) {
            for (int i = 0; i < h; i++) {
                tab[i][j] = tab[i][j-1];
            }
        }
        for (int i = 0; i < h; i++) tab[i][0] = 0;
        w++;
    }

    pair<ll, ll> getHash() {
        ll h[2];
        for (int k = 0; k < 2; k++) {
            h[k] = 0;
            for (int i = 0; i < 5; i++) {
                for (int j = 0; j < 10; j++) {
                    h[k] *= 2;
                    h[k] += tab[5*k + i][j];
                }
            }
        }

        return {h[0], h[1]};
    }

    void print() {
        for (int i = 0; i < h; i++) {
            for (int j = 0; j < w; j++) {
                cout << (tab[i][j] ? 'O' : ' ');
            }
            cout << "\n";
        }
        cout << "\n";
    }
};

vector<Tab> ps[11];

void rec(int k) {
    if (k == 1) {
        Tab t;
        t.w = 1;
        t.h = 1;
        t.tab[0][0] = 1;
        ps[k].push_back(t);
        return;
    }

    set<pair<ll, ll>> hashes;
    rec(k-1);

    for (auto &t : ps[k-1]) {
        for (int i = -1; i <= t.h; i++) {
            for (int j = -1; j <= t.w; j++) {
                if (i >= 0 && i < 10 && j >= 0 && j < 10 && t.tab[i][j] == 1) continue;
                bool any = false;
                for (int a = 0; a < 4; a++) {
                    int ai = i + dx[a];
                    int aj = j + dy[a];
                    if (ai < 0 || ai >= 10 || aj < 0 || aj >= 10 || t.tab[ai][aj] == 0) continue;
                    any = true;
                    break;
                }
                if (any) {
                    Tab nt(t);
                    if (i == -1) {
                        i++;
                        nt.shiftDown();
                    } else if (i == t.h) {
                        nt.h++;
                    }
                    if (j == -1) {
                        j++;
                        nt.shiftRight();
                    } else if (j == t.w) {
                        nt.w++;
                    }
                    nt.tab[i][j] = 1;
                    pair<ll, ll> h = nt.getHash();
                    if (hashes.find(h) == hashes.end()) {
                        ps[k].push_back(nt);
                        hashes.insert(h);
                    }
                }
            }
        }
    }
}

int cx[51][51];

void solve() {
    int n, m;
    cin >> n >> m;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cin >> cx[i][j];
        }
    }
    int res = 0;
    for (auto &p : ps[m]) {
        vector<pair<int, int>> pos;
        for (int i = 0; i < p.h; i++) {
            for (int j = 0; j < p.w; j++) {
                if (p.tab[i][j]) pos.push_back({i, j});
            }
        }
        for (int i = 0; i < n - p.h + 1; i++) {
            for (int j = 0; j < n - p.w + 1; j++) {
                int cand = 0;
                for (int k = 0; k < m; k++) {
                    int x = i + pos[k].first;
                    int y = j + pos[k].second;
                    cand += cx[x][y];
                }
                res = max(res, cand);
            }
        }
    }
    cout << res << "\n";
}

int main() {
	FASTIO;
    rec(10);
    solve();
	return 0;
}