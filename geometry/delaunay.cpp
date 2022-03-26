/* Builds a delaunay triangulation of a given set of points. This templates requires
 * the struct point defined in point_integer.cpp
 *
 * p - (input) vector of points for which the convex hull will be found.
 * adj - the delaunay triangulation of `p` as an adjacency list in which every vertex
 *       has its original index in `p`
 */

struct quadEdge {
    point o;
    quadEdge *rot, *nxt;
    bool used;
    quadEdge(point o = point(INF, INF))
        : o(o), rot(nullptr), nxt(nullptr), used(false) {}
    quadEdge* rev() const {
        return rot->rot;
    }
    quadEdge* lnext() const {
        return rot->rev()->nxt->rot;
    }
    quadEdge* prev() const {
        return rot->nxt->rot;
    }
    point dest() const {
        return rev()->o;
    }
};
 
quadEdge* makeEdge(point from, point to) {
    vector<quadEdge*> e(4);
    e[0] = new quadEdge(from);
    e[1] = new quadEdge(to);
    e[2] = new quadEdge; e[3] = new quadEdge;
    tie(e[0]->rot, e[1]->rot, e[2]->rot, e[3]->rot) = {e[2], e[3], e[1], e[0]};
    tie(e[0]->nxt, e[1]->nxt, e[2]->nxt, e[3]->nxt) = {e[0], e[1], e[3], e[2]};
    return e[0];
}
 
void splice(quadEdge* a, quadEdge* b) {
    swap(a->nxt->rot->nxt, b->nxt->rot->nxt);
    swap(a->nxt, b->nxt);
}
 
void deleteEdge(quadEdge* &e, quadEdge* ne) {
    splice(e, e->prev());
    splice(e->rev(), e->rev()->prev());
    delete e->rev()->rot;
    delete e->rev();
    delete e->rot;
    delete e;
    e = ne;
}
 
quadEdge* connect(quadEdge* a, quadEdge* b) {
    quadEdge* e = makeEdge(a->dest(), b->o);
    splice(e, a->lnext());
    splice(e->rev(), b);
    return e;
}
 
__int128 det3(point a, point b, point c) {
    vector<__int128> len = {a.norm2(), b.norm2(), c.norm2()};
    return a.x * (b.y * len[2] - c.y * len[1])
            - a.y * (b.x * len[2] - c.x * len[1])
            + len[0] * (b ^ c);
}
 
bool inCircle(point a, point b, point c, point d) {
    __int128 det = -det3(b, c, d);
    det += det3(a, c, d);
    det -= det3(a, b, d); 
    det += det3(a, b, c);
    return det > 0;
}
 
pair<quadEdge*, quadEdge*> buildTr(int l, int r, vector<point>& p) {
    if(r - l <= 3) {
        quadEdge* a = makeEdge(p[l], p[l + 1]), *b = makeEdge(p[l + 1], p[r - 1]);
        if(r - l == 2) return mp(a, a->rev());
        splice(a->rev(), b);
        ll sg = area2(p[l], p[l + 1], p[l + 2]);
        quadEdge* c = sg ? connect(b, a) : 0;
        if(sg >= 0) return mp(a, b->rev());
        else return mp(c->rev(), c);
    }
    int m = (l + r) >> 1;
    quadEdge *ldo, *ldi, *rdo, *rdi;
    tie(ldo, ldi) = buildTr(l, m, p);
    tie(rdi, rdo) = buildTr(m, r, p);
    while(1) {
        if(left(rdi->o, ldi->o, ldi->dest()))
            ldi = ldi->lnext();
        else if(right(ldi->o, rdi->o, rdi->dest()))
            rdi = rdi->rev()->nxt;
        else break;
    }
    quadEdge* basel = connect(rdi->rev(), ldi);
    auto valid = [&](quadEdge* e) {
        return right(e->dest(), basel->o, basel->dest());
    };
    if(ldi->o == ldo->o) ldo = basel->rev();
    if(rdi->o == rdo->o) rdo = basel;
    while(1) {
        quadEdge *lcand = basel->rev()->nxt;
        if(valid(lcand)) {
            while(inCircle(basel->dest(), basel->o, lcand->dest(), lcand->nxt->dest()))
                deleteEdge(lcand, lcand->nxt);
        }
        quadEdge *rcand = basel->prev();
        if(valid(rcand)) {
            while(inCircle(basel->dest(), basel->o, rcand->dest(), rcand->prev()->dest()))
                deleteEdge(rcand, rcand->prev());
        }
        if(!valid(lcand) && !valid(rcand))
            break;
        if(!valid(lcand) 
            || (valid(rcand) && inCircle(lcand->dest(), lcand->o, rcand->o, rcand->dest())))
            basel = connect(rcand, basel->rev());
        else basel = connect(basel->rev(), lcand->rev());
    }
    return mp(ldo, rdo);
}

void delaunay(vector<point> p, vector<vector<int>>& adj) {
    vector<point> temp = p;
    map<point, int> m;
    for(int i = 0; i < sz(p); i++) m[p[i]] = i;
    sort(p.begin(), p.end());
    adj.resize(sz(p));
    auto add_edge = [&](point a, point b) {
            adj[m[a]].pb(m[b]);
    };
    bool col = 1;
    for(int i = 2; i < sz(p); i++) col &= collinear(p[0], p[1], p[i]);
    if(col) {
        for(int i = 0; i + 1 < sz(p); i++) {
            add_edge(p[i], p[i + 1]);
            add_edge(p[i + 1], p[i]);
        }
    } else {
        quadEdge* e = buildTr(0, sz(p), p).first;
        vector<quadEdge*> edges = {e};
        for(int i = 0; i < sz(edges); e = edges[i++]) {
            for(quadEdge* at = e; !at->used; at = at->nxt) {
                at->used = 1;
                add_edge(at->o, at->rev()->o);
                edges.pb(at->rev());
            }
        }
    }
}