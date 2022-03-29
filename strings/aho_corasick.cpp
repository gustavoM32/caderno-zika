/* Computes the Aho-Corasick automata of a set of strings. Every node in the automata
 * corresponds to a string which its longest suffix that is a prefix of a string in the
 * set corresponds to the path of the inherent trie to that node. 
 *
 * START - first character of the alphabet 
 * ALPH - alphabet size
 *
 * Complexity: O(n * ALPH), where n is the sum of the lengths of all the strings in the set
 */

const int ALPH = 26;
const char START = 'a';
struct node {
	array<int, ALPH> go;
	int exit; 
	int link;
	node() : exit(0), link(0) {
		for(int i = 0; i < ALPH; i++)
			go[i] = 0;		
	}
};
struct ahoCorasick {
	vector<node> t;
	ahoCorasick() : t(1) {}
	inline void add(const string& s) {
		int v = 0;
		for(char c: s) {
			int cur = c - START;
			if(!t[v].go[cur]) {
				t[v].go[cur] = sz(t);
				t.pb({});
			}
			v = t[v].go[cur];
		}
		t[v].exit++;
	}
	inline void build() {
		queue<int> q;
		q.push(0);
		while(!q.empty()) {
			int v = q.front();
			q.pop();
			int l = t[v].link;
			t[v].exit += t[l].exit;
			for(int i = 0; i < ALPH; i++) {
				if(t[v].go[i]) {
					t[t[v].go[i]].link = v ? t[l].go[i] : 0;
					q.push(t[v].go[i]);
				} else t[v].go[i] = t[l].go[i];
			}
		}
	}
};