struct Median {
    multiset<ll, greater<ll>> sm;
    multiset<ll> gt;
    int size = 0;

    void balance() {
        if (sm.size() < (size + 1) / 2) {
            ll x = *(gt.begin());
            gt.erase(gt.begin());
            sm.insert(x);
        } else if (sm.size() > (size + 1) / 2) {
            ll x = *(sm.begin());
            sm.erase(sm.begin());
            gt.insert(x);
        }
    }

    void insert(ll v) {
        if (size == 0 || v <= median_first()) sm.insert(v);
        else gt.insert(v);
        size++;
        balance();
    }

    bool remove(ll v) {
        if (size == 0) return false;
        if (v <= median_first()) {
            auto it = sm.find(v);
            if (it == sm.end()) return false;
            sm.erase(it);
        } else {
            auto it = gt.find(v);
            if (it == gt.end()) return false;
            gt.erase(it);
        }
        size--;
        balance();
        return true;
    }

    ll median_first() {
        if (size % 2 == 1) return *(sm.begin());
        return *(sm.begin());
    }

    double median() {
        if (size % 2 == 1) return *(sm.begin());
        return (*(sm.begin()) + *(gt.begin())) / 2.0;
    }

    string median_string() {
        if (size % 2 == 1) return to_string(*(sm.begin()));
        ll a = *(sm.begin()) + *(gt.begin());
        string s;
        if (a < 0) s += '-';
        s += to_string(abs(a) / 2);
        if (abs(a) % 2 == 1) s += ".5";
        return s;
    }
};
