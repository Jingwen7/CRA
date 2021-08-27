#ifndef SAMPLE_H_
#define SAMPLE_H_

#include "cluster.h"
#include "index.h"
#include "graph.h"
#include "genome.h"

class sample 
{
public:
	vector<cluster> dense_clusts;
	vector<cluster> sparse_clusts;
	vector<uint32_t> breakpoints; // TODO(Jingwen) change vector to set
	vector<uint32_t> left_bps;
	vector<uint32_t> right_bps;
	vector<uint32_t> hyperbreakpoints;
	vector<uint32_t> left_hyperbps;
	vector<uint32_t> right_hyperbps;
	string *readname;
	uint32_t idx;
	uint32_t s;
	uint32_t e; // s, e are referring to the positions in superGraph.intvs
	// bool relative_strand;
	Genome *genome;

	sample () {};
	~sample () {};

	void init (uint32_t i, string *rn, Genome *g);

	void clusterBreakpoints (const fragopt_t &fopts);

	void process (vector<sample> &samples, const uint32_t a, const uint32_t b, const idx_t &mi_s, const idx_t &mi_l, 
		const fragopt_t &fopts, const idxopt_t &siopt, const idxopt_t &liopt, bool self, bool dense);

	void dump (string * readname, const fragopt_t &fopts);

	void substractClusters (bool dense);
	
	void unifyIntv(const fragopt_t &fopts, graph &superGraph, bool dense);

	void clusterBreakpointswithinInitalBreakpoints (const fragopt_t &fopts);


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

	void across_process (sample &sample_a, sample &sample_b, const uint32_t a, const uint32_t b, const idx_t &mi_s, const idx_t &mi_l, 
				const fragopt_t &fopts, const idxopt_t &siopt, const idxopt_t &liopt, bool self, bool dense);

	void decideStrand ();

	void dump (string * readname_i, string * readname_j, const fragopt_t &fopts, uint32_t len, bool dense);		

	void unifyIntv(const fragopt_t &fopts, graph &superGraph, bool dense);

	void trimclusters(bool dense);
};

void Cluster_helper (uint32_t s, uint32_t e, const vector<uint32_t> &bps, vector<uint32_t> &trimInfo, vector<uint32_t> &right, vector<uint32_t> &left, const fragopt_t &fopts);

void dumpGraph (vector<sample> &samples, graph &G, const fragopt_t &fopts, bool dense);

void FindBasicRepeatUnit(vector<sample> &samples, vector<vector<acrosample>> &acrosamples, Genome &genome, const idx_t &mi_s, const idx_t &mi_l, 
						const idxopt_t &siopt, const idxopt_t &liopt, const fragopt_t &fopts);

void FindsmallerRepeatUnit (vector<sample> &samples, vector<vector<acrosample>> &acrosamples, Genome &genome, const idx_t &mi_s, const idx_t &mi_l, 
						const idxopt_t &siopt, const idxopt_t &liopt, const fragopt_t &fopts);

void graphConstruction (graph &G, uint32_t n, vector<sample> &samples, vector<vector<acrosample>> &acrosamples, const fragopt_t &fopts, bool dense);

#endif