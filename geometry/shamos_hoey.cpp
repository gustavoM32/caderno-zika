
/* Shamos-Hoey algorithm for checking wether a collection of segments have an intersection.
 * This template depends on point_integer.cpp  and line_integer.cpp.
 *
 * seg - (input) collection of segments.
 *
 * Time complexity: O(n logn)
 */

bool shamos_hoey(vector<line> seg) {
	// create sweep line events {x, type, seg_id}
	vector<array<ll, 3> > ev;
	for(int i = 0; i < sz(seg); i++) {
	if(seg[i].q < seg[i].p) swap(seg[i].p, seg[i].q);
		ev.pb({seg[i].p.x, 0, i});
		ev.pb({seg[i].q.x, 1, i});
	}
	sort(ev.begin(), ev.end());
	set<line> s;
	for(auto e: ev) {
		line at = seg[e[2]];
		if(!e[1]) {
			auto nxt = s.lower_bound(at);
			if((nxt != s.end() && checkInter(*nxt, at))
				|| (nxt != s.begin() && checkInter(*(--nxt), at)))
					return 1;
			s.insert(at);
		} else {
			auto nxt = s.upper_bound(at), cur = nxt, prev = --cur;
			if(nxt != s.end() && prev != s.begin() 
				&& checkInter(*nxt, *(--prev))) return 1;
			s.erase(cur);
		}
	}
	return 0;
}