/* Level 1 of Segment tree beats. Can do range updates that assign to
 * every element in the range the minimum between its current value
 * and a given value x
 *
 * A - (input) initial array 
 *
 * Complexity: O(logn) per query
 *			   amortized O(logn) per update
 */

struct node{
	ll maxi_count, maxi, second, sum;
	bool leaf, lazy;
} t[N*4];

int A[N];

void build(int node) {
	t[node].maxi = max(t[node<<1].maxi, t[node<<1|1].maxi);
	t[node].sum = t[node<<1].sum + t[node<<1|1].sum;
	if(t[node<<1].maxi == t[node<<1|1].maxi) {
		t[node].maxi_count = t[node<<1].maxi_count + t[node<<1|1].maxi_count;
		t[node].second = max(t[node<<1].second, t[node<<1|1].second);
	} else if(t[node<<1].maxi > t[node<<1|1].maxi) {
		t[node].maxi_count = t[node<<1].maxi_count;
		t[node].second = max(t[node<<1|1].maxi, t[node<<1].second);
	} else {
		t[node].maxi_count = t[node<<1|1].maxi_count;
		t[node].second = max(t[node<<1].maxi, t[node<<1|1].second);
	}
}

void init(int l, int r, int node=1) {
	t[node].lazy = 0;
	t[node].leaf = 0;
	if(l == r) {
		t[node].maxi_count = 1;
		t[node].sum = t[node].maxi = A[l];
		t[node].second = -1;
		t[node].leaf = 1;
		return;
	}
	int m = (l + r) >> 1;
	init(l, m, node << 1);
	init(m + 1, r, node << 1 | 1);
	build(node);
}

void putTag(int node, ll x){
	t[node].lazy =  1;
	t[node].sum -= (t[node].maxi - x) * t[node].maxi_count;
	t[node].maxi = x;
}

void propagate(int node) {
	if(!t[node].leaf){
		if(t[node << 1].maxi > t[node].maxi)
			putTag(node << 1, t[node].maxi);
		if(t[node << 1 | 1].maxi > t[node].maxi)
			putTag(node <<1 | 1, t[node].maxi);
	}
	t[node].lazy = 0;
}

// Queries maximum in closed range [ll, rr]
int queMax(int ll, int rr, int l=0, int r=n-1, int node=1) {
	if(ll <= l && r <= rr) return t[node].maxi;
	if(rr < l || r < ll) return -1;
	if(t[node].lazy) propagate(node);
	int m = (l + r) >> 1;
	return max(queMax(ll, rr, l, m, node << 1), queMax(ll, rr, m+1, r, node << 1 | 1));
}


// Queries sum of closed range [ll, rr]
ll queSum(int ll, int rr, int l=0, int r=n-1, int node=1) {
	if(ll <= l && r <= rr) return t[node].sum;
	if(rr < l || r < ll) return 0;
	if(t[node].lazy) propagate(node);
	int m = (l + r) >> 1;
	return queSum(ll, rr, l, m, node << 1) + queSum(ll, rr, m + 1, r, node << 1 | 1);
}


// Updates to the minimum between x and current value in the closed range [ll, rr]
void upd(int ll, int rr, long long int x, int l=0, int r=n-1, int node=1) {
	if(rr < l || r < ll || t[node].maxi <= x)
		return;
	if(ll <= l && r <= rr && t[node].second < x) {
		putTag(node, x);
		return;
	}
	if(t[node].lazy) propagate(node);
	int m = (l + r) >> 1;
	upd(ll, rr, x, l, m, node << 1);
	upd(ll, rr, x, m+1, r, node << 1 | 1);
	build(node);
}
