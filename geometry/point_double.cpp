/* Basic structure of point and operations related with it. This template works
 * with double coordinates.
 *
 * All operations' time complexity are O(1)
 */

typedef long double ld;
const ld EPS = 1e-9;

struct point { 
    ld x, y;
    point(ld x, ld y) : x(x), y(y) {}
    point() {}
    ld norm2() { 
        return *this * *this;
    }
    ld norm() {
        return sqrt(norm2());
    }
    bool operator==(const point& other) const {
        return abs(x - other.x) < EPS && abs(y - other.y) < EPS;
    }
    point operator+(const point& other) const {
        return point(x + other.x, y + other.y);
    }
    point operator-(const point& other) const {
        return point(x - other.x, y - other.y);
    }
    point operator*(ld t) const {
        return point(x * t, y * t);
    }
    point operator/(ld t) const {
        return point(x / t, y / t);
    }
    ld operator*(const point& other) const {
        return x*other.x + y*other.y;
    }
    ld operator^(const point& other) const { // cross product
        return x*other.y - y*other.x;
    }
    bool operator<(const point& other) const { // for sweep line
        return x < other.x - EPS || (abs(x - other.x) < EPS && y < other.y - EPS);
    }
    point rotate(point r) {
        return point(*this ^ r, *this * r);
    }
    point rotate(ld ang) {
        return rotate(point(sin(ang), cos(ang)));
    }
    ld angle(point& other) { // only works for angles in the range [0, PI]
        ld cos_val = min(1.0L, max(-1.0L, *this * other / (norm() * other.norm())));
        return acos(cos_val);
    }
};
point ccw90(1, 0);
point cw90(-1, 0);

ld dist2(point p, point q) { // squared distance
    return (p - q).norm2();
}

ld dist(point p, point q) {
    return sqrt(dist2(p, q));
}

ld area2(point a, point b, point c) { // two times signed area of triangle abc
	return (b - a) ^ (c - a);
}

bool left(point a, point b, point c) {
	return area2(a, b, c) > EPS;
}

bool right(point a, point b, point c) {
	return area2(a, b, c) < -EPS;
}

bool collinear(point a, point b, point c) {
	return abs(area2(a, b, c)) < EPS;
}

// Returns 0 if vectors a and b are not parallel.
// If they are parallel, returns 1 if they have the same direction 
// and returns -1 otherwise
int paral(point a, point b) { 
    if((a ^ b) != 0) return 0;
    if((a.x > EPS) == (b.x > EPS) && (a.y > EPS) == (b.y > EPS))
        return 1;
    return -1;
}