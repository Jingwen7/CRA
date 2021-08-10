#include "index.h"
#include "option.h"
#include "rfpriv.h"
#include "hit.h"
#include "cluster.h"
#include "breakpoint.h"
#include "sample.h"

void sample::process (const uint32_t a, const uint32_t b, const idx_t &mi_s, const idx_t &mi_l, const fragopt_t &fopts, 
					const idxopt_t &siopt, const idxopt_t &liopt, bool self)
{
	// get diagonal cluster matches
	vector<hit> dense_fhits, dense_rhits;
	rf_hit(a, b, mi_s, dense_fhits, dense_rhits, fopts, siopt, 1); // find rhits also
	cerr << "finish rf_hit for dense kmer!" << endl;
	// clusters dense_clusts;
	cleanDiag(dense_fhits, dense_clusts, fopts, siopt, 0);
	cleanDiag(dense_rhits, dense_clusts, fopts, siopt, 1);
	dense_fhits.clear();
	dense_rhits.clear();
	cerr << "finish cleanDiag and store diagonal clusters for sparse kmer!" << endl;

	vector<hit> sparse_fhits, sparse_rhits;
	rf_hit(a, b, mi_l, sparse_fhits, sparse_rhits, fopts, liopt, 1); // find rhits also
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
	trimOnY(dense_clusts, sparse_clusts, bps_Y, fopts);

	// trim the breakpoints on x-axis;
	trimOnX(dense_clusts, sparse_clusts, bps_Y, breakpoints, fopts);
	bps_Y.clear();
}

