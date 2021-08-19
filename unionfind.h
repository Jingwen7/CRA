#ifndef UNIONFIND_H_
#define UNIONFIND_H_

#include <numeric>

using std::iota;

class UF    
{
public:
    uint32_t *id, cnt, *sz;

	// Create an empty union find data structure with N isolated sets.
    UF (uint32_t N) {
        // two dummy variables: id[0]: strand 0; id[1]: strand 1
        cnt = N + 2; // # of disjoint sets
		id = new uint32_t[N + 2];
		sz = new uint32_t[N + 2];

        for(int i = 0; i < N; i++) { 
            id[i] = i;
	    	sz[i] = 1;
        }
    }
    ~UF () {
	delete [] id;
	delete [] sz;
    }

	// Return the id of component corresponding to object p.
    uint32_t find(uint32_t p)	{
        // find the root
        uint32_t root = p;
        while (root != id[root])
            root = id[root];

        // path compression
        while (p != root) {
            uint32_t temp = id[p];
            id[p] = root;
            p = temp;
        }
        return root;
    }

	// Replace sets containing x and y with their union.
    void merge(uint32_t x, uint32_t y) {
        uint32_t i = find(x);
        uint32_t j = find(y);
        if (i == j) return;
		
		// make smaller root point to larger one
        if (sz[i] < sz[j]) { 
			id[i] = j; 
			sz[j] += sz[i]; 
		} 
		else { 
			id[j] = i; 
			sz[i] += sz[j]; 
		}
        cnt--;
    }

	// Are objects x and y in the same set?
    bool IsConnected(uint32_t x, uint32_t y) {
        return find(x) == find(y);
    }

	// Return the number of disjoint sets.
    uint32_t count() {
        return cnt;
    }
};

#endif
