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
	uint32_t s_left;
	uint32_t s_right;
	uint32_t e_left;
	uint32_t e_right;
	uint32_t sample_idx;

	interval (uint32_t S_left, uint32_t S_right, uint32_t E_left, uint32_t E_right, uint32_t Idx) : s_left(S_left), s_right(S_right), e_left(E_left), e_right(E_right), sample_idx(Idx) {};
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
	// UF * uf;
	uint32_t nofComs;

	graph () {};
	~graph () {
		// delete uf;
	};

	void init ();
	void findConnetedComponents (bool dense);
	void insertInvt (vector<uint32_t> &left_bps, vector<uint32_t> &right_bps, uint32_t idx, uint32_t &s, uint32_t &e);

};
#endif