/* Basic structure of line defined by two 2D points which the line goes through.
 * Some of the functions assume that the line is in fact a segment with the two 
 * points as its endpoints. Those operations are preceded by a comment stating 
 * that this is the case.
 *
 * This template depends on point_double.cpp and hence the coordinates don't need
 * to be integers.
 *
 * All operations' time complexity are O(1)
 */

struct line {
	point p, q;
	line(point p, point q) : p(p), pq(q) {}
	line() {}
	bool has(point r) {
        return paral((r - p), (q - p));
    }
	bool operator==(const line& other) const { // assumes that direction does not matter
        return has(other.p) && has(other.q);
    }
    bool isVert() {
        return abs(p.x - q.x) <= EPS;
    }
    point proj(point r) {
        point q_vec = q - p, r_vec = r - p;
        return p + q_vec * (q_vec * r_vec / q_vec.norm2());
    }
    ld dist(pt r) {
        return (r - proj(r)).norm();
    }
    // the following operations are for segments only
    bool segHas(point r) {
        return collinear(p, q, r)
            && (min(p.x, q.x) < r.x + EPS && r.x < max(p.x, q.x) + EPS)
            && (min(p.y, q.y) < r.y + EPS && r.y < max(p.y, q.y) + EPS);
    }
	line rotate(auto a) { // rotates segment pivoted in p
        return line(p, p + (q - p).rotate(a));
    }
};

bool paraline(line a, line b) {
    return paral(a.q - a.p, b.q - b.p);
}

point intersect(line a, line b) {
    if(paraline(a, b)) return point(INF, INF);
    point v_a = (a.q - a.p), v_b = (b.q - b.p);
    ll c_a = v_a ^ a.p, c_b = v_b ^ b.p;
    return (v_b*c_a - v_a*c_b) / (v_a ^ v_b);
}

// the following functions are for segments only
bool checkInter(line a, line b) {
    if(a.segHas(b.p) || a.segHas(b.q) || b.segHas(a.p) || b.segHas(a.q))
        return 1;
    return left(a.p, a.q, b.p) != left(a.p, a.q, b.q)
            && left(b.p, b.q, a.p) != left(b.p, b.q, a.q);
}