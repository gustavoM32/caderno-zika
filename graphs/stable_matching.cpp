
/* Solves the rural hospital version of the stable matching problem.
 * This is, you have a set of doctors and hospitals and you want to
 * assign at most cap[i] doctors to the i-th hospital and each doctor
 * to at most one hospital. Also, each doctor has a list of candidate
 * hospitals, sorted in decreasing order of preference, and each hospital
 * has a list of candidate doctors, sorted in decreasing order of 
 * preference too. Find an assigment of doctors to hospitals such that
 * there is no doctor D assigned to hospital H in which *both* D and H
 * would rather be matched with other candidates (consider the 'empty'
 * assignment as the least prefered one).
 * Remark: any stable matching has the same set of doctors assigned to
 * any hospital, and the same capactity used in every hospital.
 * 
 * For the less constrained version of the problem, just set all the
 * capacities to 1.
 *
 * n - (input) number of doctors.
 * m - (input) number of hospitals.
 * cap - (input) array of hospitals' capacities.
 * doc - (input) list of hospital candidates of each doctor sorted by preference.
 * hos - (input) list of doctor candidates of each hospital sorted by preference.
 * assign - assign[i] is the hospital assigned to doctor i or -1 if it hasn't any
 *          hospital assigned.
 *
 * Complexity: O(L), where L is the total size of the lists of candidates.
 */

int n, m;
vector<int> doc[N], hos[N];
int pos[N], cap[N], last[N];
unordered_map<int, int> prior[N];
 
vector<int> galeShapley() {
	for(int i = 0; i < m; i++) {
		last[i] = sz(hos[i]) - 1;
		for(int j = 0; j < sz(hos[i]); j++)
			prior[i][hos[i][j]] = j;
	}
	vector<int> assign(n, -1);
	queue<int> Q;
	auto tryPush = [&] (int i) {
		if(pos[i] < sz(doc[i]))
			Q.push(i);
	};
	for(int i = 0; i < n; i++) 
        tryPush(i);
	while(!Q.empty()) {
		int u = Q.front();
		Q.pop();
		int cur = doc[u][pos[u]++];
		if(cap[cur]) {
			cap[cur]--;
			assign[u] = cur;
		} else {
			if(prior[cur][u] <= last[cur]) {
				assign[u] = cur;
				while(assign[hos[cur][last[cur]]] != cur)
					last[cur]--;
				assign[hos[cur][last[cur]]] = -1;
				tryPush(hos[cur][last[cur]]);
			} else tryPush(u);
		}
	}
	return assign;
}