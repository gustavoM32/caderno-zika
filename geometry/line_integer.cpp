/* Basic structure of line defined by two 2D points which the line goes through.
 * Some of the functions assume that the line is in fact a segment with the two 
 * points as its endpoints. Those operations are preceded by a comment stating 
 * that this is the case.
 *
 * This template depends on point_integer.cpp and works only with integers.
 *
 * All operations' time complexity are O(1)
 */

struct line {
    point p, q;
    line(point p, point q) : p(p), q(q) {}
    line() {}
    bool has(const point& r) const {
        return paral((r - p), (q - p));
    }
	bool operator==(const line& other) const { // assumes that direction does not matter
        return has(other.p) && has(other.q);
    }
    bool isVert() {
        return p.x == q.x;
    }
    // the following operations are for segments only
    bool segHas(point r) {
        return collinear(p, q, r)
            && (min(p.x, q.x) <= r.x && r.x <= max(p.x, q.x))
            && (min(p.y, q.y) <= r.y && r.y <= max(p.y, q.y));
    }
    line rotate(point r) { // rotates segment pivoted in p
        return line(p, p + (q - p).rotate(r));
    }
    bool operator<(const line& other) const { // for Shamos-Hoey
        if(p == other.p) return left(p, q, other.q);
        if(!isVert() && (other.isVert() || p.x < other.p.x))
            return left(p, q, other.p);
        return left(p, other.q, other.p);
    }
};

int paraline(line a, line b) {
    return paral(a.q - a.p, b.q - b.p);
}

// the following functions are for segments only
bool checkInter(line a, line b) {
    if(a.segHas(b.p) || a.segHas(b.q) || b.segHas(a.p) || b.segHas(a.q))
        return 1;
    return left(a.p, a.q, b.p) != left(a.p, a.q, b.q)
            && left(b.p, b.q, a.p) != left(b.p, b.q, a.q);
}