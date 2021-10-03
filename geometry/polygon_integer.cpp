/* Basic structure of polygon.
 *
 * This template depends on point_integer.cpp since it only works for integer coordinates.
 *
 * All operations' time complexity are O(1) unless stated otherwise.
 */

struct polygon {
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
    ll area2() { // O(n)
        return abs(orientation());
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
    ll calipers() { // O(n). The polygon must be convex and in ccw order.
        ll ans = 0;
        for(int i = 0, j = 1; i < n; i++) {
            point vec_i = p[(i+1)%n] - p[i];
            while((vec_i ^ (p[(j+1)%n] - p[j])) > 0) 
                j = (j + 1) % n;
            ans = max(ans, dist2(p[i], p[j])); // Example with polygon diameter squared
        }
        return ans;
    }
    int extreme(const function<bool(point, point)> &cmp) {
        auto isExtreme = [&](int i, bool& curDir) -> bool {
            curDir = cmp(p[(i + 1) % n], p[i]);
            return !cmp(p[(i + n - 1) % n], p[i]) && !curDir;
        };
        bool lastDir, curDir;
        if(isExtreme(0, lastDir)) return 0;
        int lo = 0, hi = n; 
        while(lo + 1 < hi) {
            int m = (lo + hi) >> 1;
            if(isExtreme(m, curDir)) return m;
            bool relDir = cmp(p[m], p[lo]);
            if((!lastDir && curDir) || (lastDir == curDir && relDir == curDir)) {
                lo = m;
                lastDir = curDir;
            } else hi = m;
        }
        return lo;
    }
    pair<int, int> tangent(point q) { // O(log n) for convex polygon in ccw orientation
        // Finds the indices of the two tangents to an external point q
        auto leftTangent = [&](point r, point s) -> bool {
            return right(q, r, s);
        };
        auto rightTangent = [&](point r, point s) -> bool {
            return left(q, r, s);
        };
        return {extreme(leftTangent), extreme(rightTangent)};
    }
    int maximize(point v) { // O(log n) for convex polygon in ccw orientation
        // Finds the extreme point in the direction of the vector
        return extreme([&](point p, point q) {return p * v > q * v;});
    }
};
