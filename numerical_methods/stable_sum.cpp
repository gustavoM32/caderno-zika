/* From Handbook of geometry for competitive programmers - Victor Lecomte 
 * This is an algorithm to make sums of positive numbers more precise
 * Complexity is O(n) amortized and the relative error of the sum is
 * 2 log_2(n) eps, down from n * eps precision of the serial sum
 * eps is the machine precision.
 */

struct stableSum {
    int cnt = 0;
    vector<double> v, pref{0};

    // add number a to the sum
    void operator+=(double a) {
        assert(a >= 0);
        int s = ++cnt;
        while (s % 2 == 0) {
            a += v.back();
            v.pop_back();
            pref.pop_back();
            s /= 2;
        }
        v.push_back(a);
        pref.push_back(pref.back() + 1);
    }

    // return the sum value
    double val() {
        return pref.back();
    }
}