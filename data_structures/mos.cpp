/* Optimization for problems with queries that can be answered offline
 *
 * Complexity: O(q * (BLOCK_SIZE) + n * (n / BLOCK_SIZE))
 *             O((q + n) * sqrt(n)) if BLOCK_SIZE ~ sqrt(n)
 */

const int BLOCK_SIZE = 700;

struct Query {
    int l, r, idx;
    bool operator<(Query other) const
    {
        return make_pair(l / BLOCK_SIZE, r) <
               make_pair(other.l / BLOCK_SIZE, other.r);
    }
};

struct Mos {
    // variables specific for the problem
    ll sum;
    vector<ll> v;

    // stores the answers to the queries in the original order
    void exec(vector<Query> &queries, vector<ll> &answers) {
        answers.resize(queries.size());
        sort(queries.begin(), queries.end());

        int cur_l = 0;
        int cur_r = -1;

        for (Query q : queries) {
            while (cur_l > q.l) {
                cur_l--;
                add(cur_l);
            }
            while (cur_r < q.r) {
                cur_r++;
                add(cur_r);
            }
            while (cur_l < q.l) {
                remove(cur_l);
                cur_l++;
            }
            while (cur_r > q.r) {
                remove(cur_r);
                cur_r--;
            }
            answers[q.idx] = get_answer(cur_l, cur_r);
        }
    }

    // functions below are specific for the problem
    Mos(vector<ll> &v) : sum(0), v(v) {}

    void add(int i) {
        sum += v[i];
    }

    void remove(int i) {
        sum -= v[i];
    }

    ll get_answer(int l, int r) {
        return sum;
    }
};
