/* Convex hull trick version for monotonic slopes and queries, this is, 
 * the inserted lines have only increasing/decreasing slopes and the 
 * queries are done only in increasing/decreasing x-coordinates.
 *
 * This template assumes non-decreasing slopes in inserted lines
 * and queries in non-decreasing x-coordinates
 *
 * Complexity: amortized O(1) for each insertion and query
 */

typedef long double ld;

struct chullTrick{
	deque<pair<ll, ll> > lines;

    ll eval(ll x, pair<ll, ll> line){
        return line.first*x + line.second;
    }

    ld inter(pair<ll, ll> a, pair<ll, ll> b){
        return ld(b.second - a.second) / (a.first - b.first);
    }

	ll que(ll x){
		while(sz(lines) >= 2 && eval(x,lines[0]) <= eval(x,lines[1]))
			lines.pop_front();
		return eval(x, lines[0]);
	}

	void insert(pair<ll, ll> nline){
		while(sz(lines) >= 2 && inter(nline,lines[sz(lines)-2]) < inter(lines.back(),lines[sz(lines)-2]) + EPS)
			lines.pop_back();
		if(sz(lines) == 1 && lines.back().first == nline.first && lines.back().second < nline.second)
			lines.pop_back();
		if(lines.empty() || nline.first != lines.back().first)
			lines.push_back(nline);
	}	
};