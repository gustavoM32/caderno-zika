# ~/.vimrc
```vim
set nu rnu sw=4 ts=4 ai ls=2
```

# template.cpp
```cpp
#include <bits/stdc++.h>
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

void solve() {

}

signed main() {
	fastio;
	solve();
	return 0;
}
```

# Makefile
```Makefile
c = g++ -Wall -std=c++17 -static -lm $< -o $*
mtime = /usr/bin/time -f '%C %Us %MKB

%: %.cpp
	$c -g

t%: %.cpp
	$c -O2
	@for i in $*.in*; do echo "\n== $$i ==" && $(mtime) ./$* < $$i; done
```
