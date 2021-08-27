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
	// uf = new UF(intvs.size());
	uint32_t c;
	for (c = 0; c < intvs.size(); ++c) 
		lens[c] = (intvs[c].e_right + intvs[c].e_left) / 2 - (intvs[c].s_right + intvs[c].s_left) / 2;
	// intvs[c].e_right - intvs[c].s_right;
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
	cerr << "clusters: " << nofComs << endl;
}

void graph::findConnetedComponents (bool dense)
{
	if (intvs.size() == 0) return;
	nofComs = 0; 
	if (dense)
		DFS (adlist_dense, colors_dense, nofComs);
	else
		DFS (adlist_sparse, colors_sparse, nofComs);
}

void graph::insertInvt (vector<uint32_t> &left_bps, vector<uint32_t> &right_bps, uint32_t idx, uint32_t &s, uint32_t &e)
{
	// insert intervals to graph
	uint32_t i;
	uint32_t o = intvs.size();
	for (i = 1; i < left_bps.size(); ++i) 
		intvs.push_back(interval(left_bps[i - 1], right_bps[i - 1], left_bps[i], right_bps[i], idx));

	s = o;
	e = intvs.size();
}

// void graph::insertInvt (const vector<uint32_t> &clusterPivots, uint32_t idx, uint32_t &s, uint32_t &e)
// {
// 	// insert intervals to graph
// 	uint32_t i;
// 	uint32_t o = intvs.size();
// 	for (i = 1; i < clusterPivots.size(); ++i) 
// 		intvs.push_back(interval(clusterPivots[i - 1], clusterPivots[i], idx));

// 	s = o;
// 	e = intvs.size();
// }


