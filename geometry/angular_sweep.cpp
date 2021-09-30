/* Version of point_integer with < operator for angular sweep. 
 * The angle range is ]-\pi, \pi] and the template considers that the angle of the point
 * (0, 0) is 0. This template only contains a the operators needed for the for applying
 * the angular sort.
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
        return y > 0 || (y == 0 && x < 0);
    }
    bool operator==(const point& other) const {
        return x == other.x && y == other.y;
    }
	bool operator<(const point& other) const { // for angular sweep
		int this_side = side(), other_side = other.side();
        if(this_side != other_side) return this_side < other_side;
        if(*this == point(0, 0)) return 0;
        if(other == point(0, 0)) return 1;
        return (*this ^ other) > 0;
    }
};
