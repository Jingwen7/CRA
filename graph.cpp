#include "graph.h"
#include "unionfind.h"
#include <vector>
#include <set>
#include <iostream>


using std::vector;
using std::set;
using std::endl;
using std::cerr;

void graph::init ()
{
	adlist_dense.resize(intvs.size());
	adlist_sparse.resize(intvs.size());
	colors_dense.resize(intvs.size(), 0);
	colors_sparse.resize(intvs.size(), 0);
	lens.resize(intvs.size(), 0);
	uf = new UF(intvs.size());
	uint32_t c;
	for (c = 0; c < intvs.size(); ++c) 
		lens[c] = intvs[c].e - intvs[c].s;
}

void DFS (const vector<vector<uint32_t>> &adlist, vector<uint32_t> &colors, uint32_t &nofComs)
{
	vector<uint32_t> stack;
	uint32_t cur, c;
	vector<bool> visited(adlist.size(), 0);

	for (c = 0; c < adlist.size(); ++c) {
		if (visited[c])
			continue;
		nofComs++;
		cerr << "start to find a new component  " << nofComs << endl;
		// start to find a new component
		stack.push_back(c);
		colors[c] = nofComs;

		while (!stack.empty()) {
			cur = stack.back();
			stack.pop_back();
			if (adlist[cur].empty()) continue;
			for (auto&it : adlist[cur]) {
				if (visited[it]) continue;
				stack.push_back(it);
				colors[it] = nofComs;
				visited[it] = 1;
			}
		}
	}	
}

void graph::findConnetedComponents ()
{
	if (intvs.size() == 0) return;
	nofComs = 0; 
	DFS (adlist_dense, colors_dense, nofComs);
	cerr << "dense clusters: " << nofComs << endl;
	DFS (adlist_sparse, colors_sparse, nofComs);
}

