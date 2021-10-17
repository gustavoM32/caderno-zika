/* Online algorithm that computes the suffix automata of a string. 
 * Every state of the suffix automata represents a set of end positions of some substrings
 * and every path on the suffix automata represents a different substring.
 * You can use the extend method to add character by character to the suffix automata
 * or use the constructor with a string to build the suffix automara for the whole
 * string at once.
 * 
 * Complexity: O(n log ALPH) time for building the suffix automata for n characters.
 *             O(n) space
 */

struct suffixAutomata {
    struct state {
        int len, link;
        map<char, int> next;
        state() : len(0), link(-1) {}
        state(int len) : len(len) {}
        state(int len, int link, map<char, int> next)
            : len(len), link(link), next(next) {}
    };
    vector<state> st;
    int last;
    suffixAutomata() : last(0) {
        st.pb({});
    }
    suffixAutomata(string &s) : last(0) {
        st.pb({});
        for(auto c: s)
            extend(c);
    }
    void extend(char c) {
        st.pb({st[last].len + 1});
        int cur = sz(st) - 1;
        int p = last;
        while(p != -1 && !st[p].next.count(c)) {
            st[p].next[c] = cur;
            p = st[p].link;
        }
        if(p == -1) st[cur].link = 0;
        else {
            int q = st[p].next[c];
            if(st[p].len + 1 == st[q].len) st[cur].link = q;
            else {
                int clone = sz(st);
                st.pb({st[p].len + 1, st[q].link, st[q].next});
                while(p != -1 && st[p].next[c] == q) {
                    st[p].next[c] = clone;
                    p = st[p].link;
                }
                st[q].link = st[cur].link = clone;
            }
        }
        last = cur;
    }
};