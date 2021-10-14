/* Computes the Suffix Array of a set of strings.
 *
 * v - (input) vector of strings for which the suffix array will be built.
 * s - concatenation of the strings in v divided by the character '$'
 * n - suffix array size.
 * sa - suffix array. sa[i] means that the suffix s[sa[i], n) is the lexicografical i-th 
 *      suffix of s
 * lcp - array of longest common prefixes. Must call computeLcp to calculate the entries
 *      of this array. lcp[i] means that the longest common prefix between s[sa[i-1], n)
 *      and s[sa[i], n) has length lcp[i]
 * 
 * Complexity: O(n logn) time for building the suffix array.
 *             O(n) time for computing the lcp array.
 *             O(n) memory.
 */

struct suffixArray {
    string s;
    vector<int> sa, lcp;
    int n;
    void cSort(int k, vector<int>& ra) {
        int maxi = max(300, n);
        vector<int> c(maxi, 0), temp_sa(n);
        for(int i = 0; i < n; i++)
            c[i + k < n ? ra[i + k] : 0]++;
        for(int i = 1; i < maxi; i++)
            c[i] += c[i - 1];
        for(int i = n - 1; i >= 0; i--)
            temp_sa[--c[sa[i] + k < n ? ra[sa[i] + k] : 0]] = sa[i];
        sa.swap(temp_sa);
    }
    suffixArray() {}
    suffixArray(vector<string>& v) {
        for(auto str: v) {
            s += str;
            s += '$';
        }
        n = sz(s);
        sa.resize(n);
        vector<int> ra(n), temp_ra(n);
        for(int i = 0; i < n; i++) {
            ra[i] = s[i];
            sa[i] = i;
        }
        for(int k = 1; k < n; k <<= 1) {
            cSort(k, ra);
            cSort(0, ra);
            int r = temp_ra[sa[0]] = 0;
            for(int i = 1; i < n; i++)
                temp_ra[sa[i]] = (ra[sa[i]] == ra[sa[i - 1]] 
                                    && ra[sa[i] + k] == ra[sa[i - 1] + k]) ? r : ++r;
            ra.swap(temp_ra);
            if(r == n - 1) break;
        }
    }
    void computeLcp() {
        vector<int> phi(n);
        int k = 0;
        phi[sa[0]] = -1;
        for(int i = 1; i < n; i++)
            phi[sa[i]] = sa[i-1];
        vector<int> plcp(n);
        for(int i = 0; i < n - 1; i++) {
            if(phi[i] == -1) {
                plcp[i] = 0;
                continue;
            }
            while(s[i + k] == s[phi[i] + k]) k++;
            plcp[i] = k;
            k = max(k - 1, 0);
        }
        lcp.resize(n);
        for(int i = 0; i < n; i++)
            lcp[i] = plcp[sa[i]];
    }
};