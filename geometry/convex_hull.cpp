/* Finds the convex hull of a given set of points. This templates requires
 * the struct point defined in point_integer.cpp or in point_double.cpp
 *
 * p - (input) vector of points for which the convex hull will be found.
 * ch - convex hull of `p` in counter-clockwise order.
 */

vector<point> convexHull(vector<point> p) {
    int n = sz(p);
	sort(p.begin(), p.end());
	vector<point> low, up;
	for(int i = 0; i < n; i++) {
		if(i && p[i] == p[i - 1]) continue;
		while(sz(up) >= 2 && !right(up[sz(up)-2], up.back(), p[i]))
			up.pop_back();
		up.pb(p[i]);
		while(sz(low) >= 2 && !left(low[sz(low)-2], low.back(), p[i]))
			low.pop_back();
		low.pb(p[i]);
	}
    vector<point> ch;
	if(sz(low) == 1) return low;
	for(int i = 0; i < sz(low) - 1; i++)
        ch.pb(low[i]);
	for(int i = sz(up) - 1; i >= 1; i--)
        ch.pb(up[i]);
    return ch;
}