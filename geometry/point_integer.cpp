/* Basic structure of point and operations related with it. This template assumes
 * integer coordinates.
 *
 * All operations' time complexity are O(1)
 */

struct point { 
	ll x, y;
	point(ll x, ll y) : x(x), y(y){}
	point(){}
	double norm2() { 
        return *this * *this;
    }
	bool operator==(const point& other) const {
        return x == other.x && y == other.y;
    }
    point operator+(const point& other) const {
        return point(x + other.x, y + other.y);
    }
    point operator-(const point& other) const {
        return point(x - other.x, y - other.y);
    }
    point operator*(ll t) const {
        return point(x * t, y * t);
    }
    point operator/(ll t) const {
        return point(x / t, y / t);
    }
    ll operator*(const point& other) const {
        return x*other.x + y*other.y;
    }
	ll operator^(const point& other) const { // cross product
        return x*other.y - y*other.x;
    }
	bool operator<(const point& other) const { // for sweep line
		return x < other.x || (x == other.x && y < other.y);
    }
	point rotate(point r) {
        return point(*this ^ r, *this * r);
    }
};
point ccw90(1, 0);
point cw90(-1, 0);

ll dist2(point p, point q) { // squared distance
    return (p - q).norm2();
}

ll area2(point a, point b, point c) { // two times signed area of triangle abc
	return (b - a) ^ (c - a);
}

bool left(point a, point b, point c) {
	return area2(a, b, c) > 0;
}

bool right(point a, point b, point c) {
	return area2(a, b, c) < 0;
}

bool collinear(point a, point b, point c) {
	return abs(area2(a, b, c)) == 0;
}

// Returns 0 if vectors a and b are not parallel.
// If they are parallel, returns 1 if they have the same direction 
// and returns -1 otherwise
int paral(point a, point b) { 
    if(a ^ b != 0) return 0;
    if((a.x > 0) == (b.x > 0) && (a.y > 0) == (b.y > 0))
        return 1;
    return -1;
}