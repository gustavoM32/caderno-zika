/* Basic structure of polygon.
 *
 * This template depends on point_double.cpp since it can work with double coordinates.
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
    ld orientation() { // O(n)
        ld acum = 0;
        for(int i = 0; i < n; i++)
            acum += p[i] ^ p[(i + 1) % n];
        return acum;
    }
    ld area() { // O(n)
        return abs(orientation()) / 2.0;
    }
    void turnCcw() { // O(n)
        if(orientation() < -EPS)
            reverse(p.begin(), p.end());
    }
    bool has(point q) { // O(log n). The polygon must be convex and in ccw order.
        if(right(p[0], p[1], q) || left(p[0], p[n-1], q)) return 0;
        int lo = 0, hi = n;
        while(lo + 1 < hi) {
            int mid = (lo + hi) >> 1;
            if(left(p[0], p[mid], q)) lo = mid;
            else hi = mid;
        }
        return !right(p[lo], p[(lo + 1) % n], q);
    }
    ld calipers() { // O(n). The polygon must be convex and in ccw order.
        ld ans = 0;
        for(int i = 0, j = 1; i < n; i++) {
            point vec_i = p[(i+1)%n] - p[i];
            while((vec_i ^ (p[(j+1)%n] - p[j])) > EPS) 
                j = (j + 1) % n;
            ans = max(ans, dist(p[i], p[j])); // Example with polygon diameter
        }
        return ans;
    }
};
