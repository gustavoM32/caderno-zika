/* Alternative version of point_integer with < operator for angular sweep. 
 * The order given by the operator is counter-clockwise starting from the ray
 * (y <= 0, x = 0). This template only contains a subset of operators and 
 * functions of the full points template.
 */

struct point { 
	ll x, y;
	point(ll x, ll y) : x(x), y(y) {}
	point() {}
    ll operator*(const point& other) const {
        return x*other.x + y*other.y;
    }
	ll operator^(const point& other) const { // cross product
        return x*other.y - y*other.x;
    }
    int side() const {
        return x < 0 || (x == 0 && y > 0);
    }
	bool operator<(const point& other) const { // for angular sweep
		int this_side = side(), other_side = other.side();
        if(this_side != other_side) return this_side < other_side;
        return *this ^ other > 0
    }
	point rotate(point r) {
        return point(*this ^ r, *this * r);
    }
};
point ccw90(1, 0);
point cw90(-1, 0);