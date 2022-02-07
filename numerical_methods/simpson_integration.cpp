/* Returns the integral of some function f in the interval [a, b].
 *
 * Complexity: O(2n)
 */
double f(double x);

double integral(double a, double b, int n) {
    n *= 2;
    double h = (b - a) / n;
    double s = f(a) + f(b);
    for (int i = 1; i < n; i++) {
        s += f(a + h*i) * (1 + (i & 1)) * 2;
    }
    return s * h / 3;
}

