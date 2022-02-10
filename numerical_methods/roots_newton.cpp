/* Find some root of a function f with derivative df.
 *
 * Complexity: The precision doubles for each iteration, in ideal conditions.
 * For roots with multiplicity greater than one, the precision increases linearly.
 */
double f(double x);
double df(double x);

double findRoot(double x0=0.0) {
    double x = x0;
    double fx = f(x);
    while (abs(fx) > EPS) {
        x -= fx / df(x);
        fx = f(x);
    }
    return x;
}
