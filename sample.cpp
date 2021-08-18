#include "index.h"
#include "option.h"
#include "rfpriv.h"
#include "hit.h"
#include "cluster.h"
#include "breakpoint.h"
#include "sample.h"

void sample::init (uint32_t i, string *rn)
{
	readname = rn;
	idx = i;
	dense_clusts.clear();
	sparse_clusts.clear();
}

void sample::process (vector<sample> &samples, const uint32_t a, const uint32_t b, const idx_t &mi_s, const idx_t &mi_l, const fragopt_t &fopts, 
					const idxopt_t &siopt, const idxopt_t &liopt, bool self)
{
	// get diagonal cluster matches
	vector<hit> dense_fhits, dense_rhits;
	rf_hit(samples, a, b, mi_s, dense_fhits, dense_rhits, fopts, siopt, 1); // find rhits also
	cerr << "finish rf_hit for dense kmer!" << endl;

	// clusters dense_clusts;
	cleanDiag(dense_fhits, dense_clusts, fopts, siopt, 0);
	cleanDiag(dense_rhits, dense_clusts, fopts, siopt, 1);
	dense_fhits.clear();
	dense_rhits.clear();
	cerr << "finish cleanDiag and store diagonal clusters for sparse kmer!" << endl;

	vector<hit> sparse_fhits, sparse_rhits;
	rf_hit(samples, a, b, mi_l, sparse_fhits, sparse_rhits, fopts, liopt, 1); // find rhits also
	cerr << "finish rf_hit for dense kmer!" << endl;

	// clusters sparse_clusts;
	cleanDiag(sparse_fhits, sparse_clusts, fopts, liopt, 0);
	cleanDiag(sparse_rhits, sparse_clusts, fopts, liopt, 1);
	sparse_fhits.clear();
	sparse_rhits.clear();
	cerr << "finish cleanDiag and store diagonal clusters for sparse kmer!" << endl;

	// flip the clusters
	flipCluster(dense_clusts, fopts, siopt, 1);
	flipCluster(sparse_clusts, fopts, liopt, 0);

	// (TODO) Jingwen: make sure the code work for inversed cluster
	// trim the breakpoints on y-axis
	vector<uint32_t> bps_Y;
	set<uint32_t> bpset;
	firstTrim(samples[a], bpset, dense_clusts, sparse_clusts, bps_Y, fopts, 0, 0);

	// trim the breakpoints on x-axis;
	secondTrim(samples[a], dense_clusts, sparse_clusts, bps_Y, breakpoints, fopts, 0, 0);
	bps_Y.clear();
}

void sample::dump (string * readname, const fragopt_t &fopts)
{
	if (fopts.debug) {
		ofstream clust("cluster.bed", ios_base::app);
		for (int m = 0; m < dense_clusts.size(); m++) {
			clust << dense_clusts[m].xStart << "\t" << dense_clusts[m].yStart << "\t" << dense_clusts[m].xEnd << "\t"
				   << dense_clusts[m].yEnd << "\t" << dense_clusts[m].xEnd - dense_clusts[m].xStart << "\t"  
				   << dense_clusts[m].strand  << "\t" << m << "\t" << "1" << "\t" << *readname << endl;				
		}
		for (int m = 0; m < sparse_clusts.size(); m++) {
			clust << sparse_clusts[m].xStart << "\t" << sparse_clusts[m].yStart << "\t" << sparse_clusts[m].xEnd << "\t"
				   << sparse_clusts[m].yEnd << "\t" << sparse_clusts[m].xEnd - sparse_clusts[m].xStart << "\t"  
				   << sparse_clusts[m].strand  << "\t" << m << "\t" << "0" << "\t" << *readname << endl;				
		}
		clust.close();
		ofstream fclust("trimlines.bed", ios_base::app);
	  	for (auto& it : breakpoints)
	    	fclust << it << "\t" << *readname << endl;
		fclust.close();			
	}
}

void sample::substractClusters ()
{
	int size = dense_clusts.size();
	dense_clusts.resize( size / 2);
	size = sparse_clusts.size();
	sparse_clusts.resize(size / 2);
}

void sample::clusterBreakpoints (const fragopt_t &fopts, graph &superGraph)
{
	vector<uint32_t> trimInfo(breakpoints.size());
	iota(trimInfo.begin(), trimInfo.end(), 0);
	vector<uint32_t> clusterPivots;

	uint32_t prev = breakpoints.back();
	int i = breakpoints.size() - 2; 
	while (i >= 0) {
		if (breakpoints[i] + fopts.clusterTrimedge >= prev)
			trimInfo[i] = trimInfo[i + 1];
		else
			trimInfo[i] = i;
		prev = breakpoints[i];
		i--;
	}

	for (i = 0; i < trimInfo.size(); ++i) {
		if (trimInfo[i] == i)
			clusterPivots.push_back(breakpoints[i]);
	}

	breakpoints = clusterPivots; // update breakpoints to cleaner version

	// insert intervals to graph
	uint32_t o = superGraph.intvs.size();
	for (i = 1; i < clusterPivots.size(); ++i) 
		superGraph.intvs.push_back(interval(clusterPivots[i - 1], clusterPivots[i], idx));

	s = o;
	e = superGraph.intvs.size();

	// o = superGraph.adlist_dense.size();
	// superGraph.adlist_dense.resize(o + clusterPivots.size());

	// o = superGraph.adlist_sparse.size();
	// superGraph.adlist_sparse.resize(o + clusterPivots.size());
}

bool overlap (uint32_t b_s, uint32_t b_e, uint32_t a_s, uint32_t a_e, double frac)
{
		uint32_t ovp = 0;

		if (b_s >= a_s and b_s < a_e) {
			ovp = min(a_e, b_e) - b_s;
		}
		else if (b_e > a_s and b_e < a_e) {
			ovp = b_e - max(a_s, b_s);
		}
		else if (b_s <= a_s and b_e > a_e) {
			ovp = a_e - a_s;
		}
		double denomA = (double) (a_e - a_s);
		double denomB = (double) (b_e - b_s);
		if ( max(ovp / denomA, ovp / denomB) >= frac) 
			return true;
		else 
			return false;
}

void addEdge(const fragopt_t &fopts, const vector<cluster> &clusts, vector<uint32_t> &breakpoints, vector<vector<uint32_t>> &adlist, uint32_t s)
{
	uint32_t i;
	vector<uint32_t> xhits, yhits;

	for (auto&it : clusts) {
		xhits.clear();
		yhits.clear();
		for (i = 1; i < breakpoints.size(); ++i) {
			if (overlap(breakpoints[i - 1], breakpoints[i], it.xStart, it.xEnd, fopts.ovpfrac))
				xhits.push_back(i - 1);
			if (overlap(breakpoints[i - 1], breakpoints[i], it.yStart, it.yEnd, fopts.ovpfrac))
				yhits.push_back(i - 1);
		}
		assert(xhits.size() == yhits.size());

		for (i = 0; i < xhits.size(); ++i) {
			adlist[xhits[i] + s].push_back(yhits[i] + s);
			adlist[yhits[i] + s].push_back(xhits[i] + s);
		}
	}
}

void sample::unifyIntv(const fragopt_t &fopts, graph &superGraph)
{
	addEdge(fopts, dense_clusts, breakpoints, superGraph.adlist_dense, s);
	addEdge(fopts, sparse_clusts, breakpoints, superGraph.adlist_sparse, s);
}

void acrosample::init (uint32_t idx, uint32_t jdx, sample * Pi, sample * Pj)
{
	i = idx;
	j = jdx;
	pi = Pi;
	pj = Pj;
	dense_clusts.clear();
	sparse_clusts.clear();
}

void acrosample::across_process (vector<sample> &samples, const uint32_t a, const uint32_t b, const idx_t &mi_s, const idx_t &mi_l, 
		const fragopt_t &fopts, const idxopt_t &siopt, const idxopt_t &liopt, bool self)
{
	cerr << "process across sample " << a << " " << b << endl;
	// get diagonal cluster matches
	vector<hit> dense_fhits, dense_rhits;
	rf_hit(samples, a, b, mi_s, dense_fhits, dense_rhits, fopts, siopt, 0); // find rhits also
	cerr << "finish rf_hit for dense kmer!" << endl;

	// clusters dense_clusts;
	cleanDiag(dense_fhits, dense_clusts, fopts, siopt, 0);
	cleanDiag(dense_rhits, dense_clusts, fopts, siopt, 1);

	dense_fhits.clear();
	dense_rhits.clear();
	cerr << "finish cleanDiag and store diagonal clusters for dense kmer!" << endl;

	vector<hit> sparse_fhits, sparse_rhits;
	rf_hit(samples, a, b, mi_l, sparse_fhits, sparse_rhits, fopts, liopt, 0); // find rhits also
	cerr << "finish rf_hit for sparse kmer!" << endl;

	// clusters sparse_clusts;
	cleanDiag(sparse_fhits, sparse_clusts, fopts, liopt, 0);
	cleanDiag(sparse_rhits, sparse_clusts, fopts, liopt, 1);
	sparse_fhits.clear();
	sparse_rhits.clear();
	cerr << "finish cleanDiag and store diagonal clusters for sparse kmer!" << endl;

	// // (TODO) Jingwen: make sure the code work for inversed cluster
	// // trim the breakpoints on y-axis
	// vector<uint32_t> bps_Y;
	// trimOnY(readname, dense_clusts, sparse_clusts, bps_Y, fopts);

	// // trim the breakpoints on x-axis;
	// trimOnX(readname, dense_clusts, sparse_clusts, bps_Y, breakpoints, fopts);
	// bps_Y.clear();
}

void acrosample::decideStrand ()
{
	uint32_t flen = 0, rlen = 0;
	uint32_t c;

	for (c = 0; c < dense_clusts.size(); ++c) {
		if (dense_clusts[c].strand == 0) 
			flen += dense_clusts[c].xEnd - dense_clusts[c].xStart;
		else
			rlen += dense_clusts[c].xEnd - dense_clusts[c].xStart;
	}

	if (flen >= rlen) 
		strand = 0;
	else
		strand = 1;
}

void acrosample::dump (string * readname_i, string * readname_j, const fragopt_t &fopts, uint32_t len)
{
	if (fopts.debug) {
		ofstream clust("cluster_acrosssample.bed", ios_base::app);
		cerr << "across sample dense_clusts.size(): " << dense_clusts.size() << endl;
		cerr << "across sample sparse_clusts.size(): " << sparse_clusts.size() << endl;

		for (int m = 0; m < dense_clusts.size(); ++m) {
			clust << dense_clusts[m].xStart << "\t" << dense_clusts[m].yStart << "\t" << dense_clusts[m].xEnd << "\t"
				   << dense_clusts[m].yEnd << "\t" << dense_clusts[m].xEnd - dense_clusts[m].xStart << "\t"  
				   << dense_clusts[m].strand  << "\t" << m << "\t" << "1" << "\t" << *readname_i << "\t" << *readname_j << endl;				
		}
		for (int m = 0; m < sparse_clusts.size(); ++m) {
			clust << sparse_clusts[m].xStart << "\t" << sparse_clusts[m].yStart << "\t" << sparse_clusts[m].xEnd << "\t"
				   << sparse_clusts[m].yEnd << "\t" << sparse_clusts[m].xEnd - sparse_clusts[m].xStart << "\t"  
				   << sparse_clusts[m].strand  << "\t" << m << "\t" << "0" << "\t" << *readname_i << "\t" << *readname_j << endl;				
		}			
		clust.close();

		ofstream fclust("trimlines_acrosssample.bed", ios_base::app);
	  	for (auto& it : pi->breakpoints)
	    	fclust << it << "\t" << *readname_i << "\t" << *readname_j << endl;	
	  	for (auto& it : pj->breakpoints)
	    	fclust << it << "\t" << *readname_j << "\t" << *readname_i << endl;
		fclust.close();			
	}
}

void addEdge_acros (const fragopt_t &fopts, const vector<cluster> &clusts, sample *pi, sample *pj, vector<vector<uint32_t>> &adlist, bool strand)
{
	uint32_t c;
	uint32_t t;
	vector<uint32_t> xhits, yhits;
	
	for (auto&it : clusts) {
		xhits.clear();
		yhits.clear();
		for (c = 1; c < pi->breakpoints.size(); ++c) {
			if (overlap(pi->breakpoints[c - 1], pi->breakpoints[c], it.xStart, it.xEnd, fopts.ovpfrac))
				xhits.push_back(c - 1);
		}

		for (c = 1; c < pj->breakpoints.size(); ++c) {
			if (overlap(pj->breakpoints[c - 1], pj->breakpoints[c], it.yStart, it.yEnd, fopts.ovpfrac))
				yhits.push_back(c - 1);
		}

		assert(xhits.size() == yhits.size());

		if (strand == 0) {
			for (c = 0; c < xhits.size(); ++c) {
				adlist[xhits[c] + pi->s].push_back(yhits[c] + pj->s);
				adlist[yhits[c] + pj->s].push_back(xhits[c] + pi->s);
			}
		}
		else {
			for (c = 0; c < xhits.size(); ++c) {
				t = xhits.size() - 1 - c;
				adlist[xhits[c] + pi->s].push_back(yhits[t] + pj->s);
				adlist[yhits[t] + pj->s].push_back(xhits[c] + pi->s);
			}			
		}
	}	
}

void acrosample::unifyIntv(const fragopt_t &fopts, graph &superGraph)
{
	cerr << "dense" << endl;
	addEdge_acros(fopts, dense_clusts, pi, pj, superGraph.adlist_dense, strand);
	cerr << "sparse" << endl;
	addEdge_acros(fopts, sparse_clusts, pi, pj, superGraph.adlist_sparse, strand);
}
















