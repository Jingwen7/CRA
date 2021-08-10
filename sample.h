#ifndef SAMPLE_H_
#define SAMPLE_H_

#include "cluster.h"

class sample 
{
public:
	vector<cluster> dense_clusts;
	vector<cluster> sparse_clusts;
	vector<uint32_t> breakpoints;

	sample () {};
	~sample () {};

	void process (const uint32_t a, const uint32_t b, const idx_t &mi_s, const idx_t &mi_l, 
		const fragopt_t &fopts, const idxopt_t &siopt, const idxopt_t &liopt, bool self);



};

#endif