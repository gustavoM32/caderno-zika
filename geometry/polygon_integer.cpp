/* Basic structure of polygon.
 *
 * This template depends on point_integer.cpp since it only works for integer coordinates.
 *
 * All operations' time complexity are O(1) unless stated otherwise.
 */

struct polygon{
	vector<point> p;
	int n;
    polygon() : n(0) {}
    polygon(vector<point> _p) {
        p = _p;
        n = sz(p);
    }
	void add(point q) {
		p.pb(q);
        n++;
	}
    // If positive, the polygon is in ccw order. It is in cw order otherwise.
    ll orientation() { // O(n)
        ll acum = 0;
        for(int i = 0; i < n; i++)
            acum += p[i] ^ p[(i + 1) % n];
        return acum;
    }
    void turnCcw() { // O(n)
        if(orientation() < 0)
            reverse(p.begin(), p.end());
    }
	bool has(point q) { // O(log n). The polygon must be convex and in ccw order
		if(right(p[0], p[1], q) || left(p[0], p[n-1], q)) return 0;
		int lo = 0, hi = n;
		while(lo + 1 < hi) {
			int mid = (lo + hi) >> 1;
			if(left(p[0], p[mid], q)) lo = mid;
			else hi = mid;
		}
		return !right(p[lo], p[(lo + 1) % n], q);
	}
};
