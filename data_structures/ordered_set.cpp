/* The ordered set data structure works just like the std::set, but it has additional
 * information for every node. Namely, it can retrieve the order for a given key and also
 * the key of a node which is in a given position.
 *
 * Methods:
 * - find_by_order(pos)
 *   returns the iterator of the node in the `pos` position.
 * - order_of_key(key)
 *   returns the order in which a node of key `key` would be in the ordered set
 * The other methods are the same of set's methods.
 *
 * Time complexity: O(logn) for insert, erase and queries
 * Space complexity: O(n)
 */

#include <ext/pb_ds/assoc_container.hpp>
#include <ext/pb_ds/tree_policy.hpp>
using namespace __gnu_pbds;

template<typename T>
using ordered_set = tree<T, null_type, less<T>, rb_tree_tag, tree_order_statistics_node_update>;