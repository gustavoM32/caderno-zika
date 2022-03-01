/* Find the minimum of a function that first strictly decreases, then has its
 * minimum and finally strictly increases
 *
 * Complexity: O(log(n))
 */

double ternary_search(double l, double r) {
	while (r - l > EPS) { // TODO set EPS
		double m1 = l + (r - l) / 3;
		double m2 = r - (r - l) / 3;
		if (f(m1) > f(m2)) l = m1; // < for maximum
		else r = m2;
	}
	return l;
}
