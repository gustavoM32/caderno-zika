/* Implicit treap with lazy propagation. It is easy to adapt to versions without
 * lazy propagation or that are not implicit.
 * For non-lazy versions, just delete the functions apply and push.
 * For non-implicit versions, use an actual key in the split instead of the implicit key.
 *
 * This example uses lazy propagation for reversing a range in implicit treaps.
 *
 * Time complexity: O(log n) for merge and split.
 * Space complexity: O(n)
 */

mt19937 rng(chrono::steady_clock::now().time_since_epoch().count());

typedef struct item* pitem;
struct item {
	int cnt, pr, val;
	bool lazy;
	pitem l, r;
	item(int val) 
        : pr(rng()), val(val), cnt(1), lazy(0), l(nullptr), r(nullptr) {}
};

int cnt(pitem t) {
	return t != nullptr ? t->cnt : 0;
}

void upd(pitem t) {
	if(t != nullptr)
		t->cnt = cnt(t->l) + 1 + cnt(t->r);
}

void apply(pitem t, bool d) {
	if(t != nullptr) {
		if(d.first) swap(t->l, t->r);
		t->lazy ^= d;
	}	
}

void push(pitem t) {
	if(t != nullptr) {
		apply(t->l, t->lazy);
		apply(t->r, t->lazy);
		t->lazy = {0, 0};
	}
}

void merge(pitem& t, pitem l, pitem r) {
	if(l == nullptr || r == nullptr)
		t = l != nullptr ? l : r;
	else if(l->pr > r->pr) {
		push(l);
		merge(l->r, l->r, r), t = l;
	} else {
		push(r);
		merge(r->l, l, r->l), t = r;
	}
	upd(t);
}

void split(pitem t, pitem& l, pitem& r, int key, int add=0) {
	if(t == nullptr) {
		l = r = nullptr; return;
	}
	push(t);
	int cur_key = add + cnt(t->l);
	if(key <= cur_key) 
		split(t->l, l, t->l, key, add), r = t;
	else
		split(t->r, t->r, r, key, add + 1 + cnt(t->l)), l = t;
	upd(t);
}

void insert(pitem& t, int pos) {
	pitem it = new item(pos);
	pitem t1, t2;
	split(t, t1, t2, pos);
	merge(t1, t1, it);
	merge(t, t1, t2);
}