/* Computes the longest (weakly) increasing subsequence in a vector.
 * Is easy to change the code below for the strictly increasing version of the problem.
 *
 * v - (input) vector for which the longest increasing subsequence will be computed.
 * aux - one of the longest increasing subsequences of v
 *
 * Complexity: O(n logn) time
 *             O(n) memory
 */

vector<int> lis(vector<int>& v) { 
    vector<int> aux;
	for(auto x: v) {
        // Change to lower_bound for strictly increasing version
		auto it = upper_bound(aux.begin(), aux.end(), x);
		if(it == aux.end()) aux.pb(x);
		else *it = x;
	}
	return aux;
}