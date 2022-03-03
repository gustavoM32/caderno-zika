/* Useful polynomials operation of competitive programming
 * From https://github.com/e-maxx-eng/e-maxx-eng-aux/blob/master/src/polynomial.cpp
 * The main operations are kept in this file.
 *
 * Complexity: varies
 */
struct Poly {
	vector<ll> p;

	Poly() {};
	
	Poly(int deg) { // for deg < 0, p(x) = 0
		p.resize(max(0, deg+1));
	}

	Poly(vector<ll> &coefs) : p(coefs) {
		cout << coefs[0] << " " << coefs[1] << endl;
	}

	int deg() { // degree is -1 for p(x) = 0
		return int(p.size()) - 1;
	}

	Poly operator*(ll x) { // O(deg())
		Poly b = *this;
		for (ll &c : b.p) c *= x;
		return b;
	}

	Poly operator*=(ll x) {
		return *this = *this * x;
	}

	Poly mul_slow(Poly &a) { // O(deg() * a.deg())
		Poly b(deg() + a.deg());
		for (int i = 0; i <= deg(); i++) {
			for (int j = 0; j <= a.deg(); j++) {
				b.p[i+j] += p[i] * a.p[j];
			}
		}
		return b;
	}

	void print() { // debug
		for (int i = deg(); i >= 0; i--) {
			if (p[i] < 0) cout << "- ";
			else cout << "+ ";
			cout << abs(p[i]);
			if (i != 0) cout << "x^" << i << " ";
		}
		cout << endl;
	}
};

