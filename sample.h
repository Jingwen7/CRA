#ifndef SAMPLE_H_
#define SAMPLE_H_

#include "cluster.h"
#include "index.h"

class sample 
{
public:
	vector<cluster> dense_clusts;
	vector<cluster> sparse_clusts;
	vector<uint32_t> breakpoints;
	string *readname;
	uint32_t idx;

	sample () {};
	~sample () {};

	void init (uint32_t i, string *rn);

	void process (vector<sample> &samples, const uint32_t a, const uint32_t b, const idx_t &mi_s, const idx_t &mi_l, 
		const fragopt_t &fopts, const idxopt_t &siopt, const idxopt_t &liopt, bool self);

	void dump (string * readname, const fragopt_t &fopts);
};


// (TODO) Jingwen: if the across sample dotplot is reverse stranded, need to reverse breakpoints of one sample
class acrosample
{
public:
	bool strand;
	uint32_t i; // assume i < j
	uint32_t j;
	sample * pi;
	sample * pj;
	vector<cluster> dense_clusts; // x-axis: sample_i; y-axis: samle-j
	vector<cluster> sparse_clusts;

	acrosample () {};
	~acrosample () {};

	void init (uint32_t idx, uint32_t jdx, sample * Pi, sample * Pj);

	void across_process (vector<sample> &samples, const uint32_t a, const uint32_t b, const idx_t &mi_s, const idx_t &mi_l, 
		const fragopt_t &fopts, const idxopt_t &siopt, const idxopt_t &liopt, bool self);
	
	void decideStrand ();

	void dump (string * readname_i, string * readname_j, const fragopt_t &fopts, uint32_t len);		
};



#endif