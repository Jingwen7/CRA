#ifndef GRAPH_H_
#define GRAPH_H_

#include <stdint.h>
#include <vector>
#include <set>
#include <iostream>
#include "unionfind.h"

using std::vector;
using std::set;
using std::endl;
using std::cerr;
using std::vector;

class interval
{
public:
	uint32_t s;
	uint32_t e;
	uint32_t sample_idx;

	interval (uint32_t S, uint32_t E, uint32_t Idx) : s(S), e(E), sample_idx(Idx) {};
	~interval () {};
};

class graph
{
public:
	vector<interval> intvs; 
	vector<vector<uint32_t>> adlist_dense; // u : vector<uint32_t> : u ---> v
	vector<vector<uint32_t>> adlist_sparse;
	vector<uint32_t> colors_dense;
	vector<uint32_t> colors_sparse;
	vector<uint32_t> lens;
	UF * uf;
	uint32_t nofComs;

	graph () {};
	~graph () {
		delete uf;
	};

	void init ();
	void findConnetedComponents ();
};
#endif