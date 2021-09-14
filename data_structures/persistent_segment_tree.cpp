/* Persistent segment tree. This example is for queries of sum in range
 * and updates of sum in position, but any query or update can be achieved
 * changing the NEUT value, and functions updNode and merge.
 * The version 0 of the persistent segTree has implicitly an array of length
 * n full of NEUT values.
 *
 * It's recommend to set n as the actual length of the array.
 * ```
 * int n; cin >> n;
 * segTree::n = n;
 * ```
 *
 * Complexity: O(logn) memory and time per query/update
 */

const int NEUT = 0;

struct segTree {
	vector<int> t = vi(1, NEUT), left = vi(1, 0), right = vi(1, 0);
	static int n;
	int newNode(int v, int l=0, int r=0) {
		t.pb(v), left.pb(l), right.pb(r);
		return sz(t) - 1;
	}
	int merge(int a, int b) {
		return a + b;
	}
	// Initializes a segTree with the values of the array A of length n
	int init(int* A, int L=0, int R=n) {
		if(L + 1 == R) return newNode(A[L]);
		int M = (L + R) >> 1;
		int l = init(A, L, M), r = init(A, M, R);
		return newNode(merge(t[l], t[r]), l , r);
	}
	int updNode(int cur_value, int upd_value) {
		return cur_value + upd_value;
	}
	int upd(int k, int pos, int v, int L=0, int R=n) {
		int nxt = newNode(t[k], left[k], right[k]);
		if(L+1 == R) t[nxt] = updNode(t[nxt], v);
		else {
			int M = (L + R) >> 1;
			int temp;
			if(pos < M) temp = upd(left[nxt], pos, v, L, M), left[nxt] = temp;
			else temp = upd(right[nxt], pos, v, M, R), right[nxt] = temp;
			t[nxt] = merge(t[left[nxt]], t[right[nxt]]);
		}
		return nxt;
	}
	int que(int k, int l, int r, int L=0, int R=n){
		if(r <= L || R <= l) return NEUT;
		if(l <= L && R <= r) return t[k];
		int M = (L + R) >> 1;
		return merge(que(left[k], l, r, L, M), que(right[k], l, r, M, R));
	}
};
int segTree::n = N;