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


















