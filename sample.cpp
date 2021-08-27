#include "index.h"
#include "option.h"
#include "rfpriv.h"
#include "hit.h"
#include "cluster.h"
#include "breakpoint.h"
#include "sample.h"
#include "unionfind.h"
#include <boost/range/adaptor/reversed.hpp> 

void sample::init (uint32_t i, string *rn, Genome *G)
{
	readname = rn;
	idx = i;
	dense_clusts.clear();
	sparse_clusts.clear();
	genome = G;
	breakpoints.clear();
	left_bps.clear();
	right_bps.clear();
}

void sample::process (vector<sample> &samples, const uint32_t a, const uint32_t b, const idx_t &mi_s, const idx_t &mi_l, const fragopt_t &fopts, 
					const idxopt_t &siopt, const idxopt_t &liopt, bool self, bool dense)
{
	// get diagonal cluster matches
	set<uint32_t> bps;
	if (dense) {
		vector<hit> dense_fhits, dense_rhits;
		rf_hit(samples[a], samples[b], a, b, mi_s, dense_fhits, dense_rhits, fopts, siopt, 1); // find rhits also
		cerr << "finish rf_hit for dense kmer!" << endl;

		// clusters dense_clusts;
		cleanDiag(dense_fhits, dense_clusts, fopts, siopt, 0);
		cleanDiag(dense_rhits, dense_clusts, fopts, siopt, 1);

		if (fopts.debug) {
			ofstream clust("single_hits.bed", ios_base::app);
			for (int m = 0; m < dense_fhits.size(); m++) {
				clust << dense_fhits[m].x << "\t" << dense_fhits[m].y << "\t" << dense_fhits[m].x + mi_s.k << "\t"
					   << dense_fhits[m].y + mi_s.k << "\t" << mi_s.k << "\t"  
					   << "0"  << "\t" << "1" << "\t" << *(samples[a].readname) << endl;				
			}
			for (int m = 0; m < dense_rhits.size(); m++) {
				clust << dense_rhits[m].x << "\t" << dense_rhits[m].y << "\t" << dense_rhits[m].x + mi_s.k << "\t"
					   << dense_rhits[m].y + mi_s.k << "\t" << mi_s.k << "\t"  
					   << "1" << "\t" << "1" << "\t" << *(samples[a].readname) << endl;				
			}
			clust.close();			
		}

		dense_fhits.clear();
		dense_rhits.clear();
		cerr << "finish cleanDiag and store diagonal clusters for dense kmer!" << endl;			

		// flip the clusters
		flipCluster(dense_clusts, fopts, siopt, 1);
		// uint32_t len = samples[a].genome->getLen(samples[a].idx);
		// dense_clusts.push_back(cluster(0, len, 0, len, 0));

		// assert(dense_clusts.size() > 0 and dense_clusts.size() % 2 == 0);

		// insert original clusters boundaries to bps
		for (auto&it : dense_clusts) {
			bps.insert(it.yStart);
			bps.insert(it.yEnd);
			bps.insert(it.xStart);
			bps.insert(it.xEnd);
		}

		// insert the original breakpoints
		if (samples[a].breakpoints.size() > 0) {
			for (auto&it : samples[a].breakpoints)
				bps.insert(it);			
		}

		samples[a].breakpoints.clear();

		if (bps.size() > 0) {
		  	for (auto&it : bps) 
		    	samples[a].breakpoints.push_back(it);			
		}


	    cerr << "sample " << a << " breakpoints.size(): " << samples[a].breakpoints.size() << endl;
	  	
	  	bps.clear();
	}
	else {
		vector<hit> sparse_fhits, sparse_rhits;
		rf_hit(samples[a], samples[b], a, b, mi_l, sparse_fhits, sparse_rhits, fopts, liopt, 1); // find rhits also
		cerr << "finish rf_hit for sparse kmer!" << endl;

		// clusters sparse_clusts;
		cleanDiag(sparse_fhits, sparse_clusts, fopts, liopt, 0);
		cleanDiag(sparse_rhits, sparse_clusts, fopts, liopt, 1);
		sparse_fhits.clear();
		sparse_rhits.clear();
		cerr << "finish cleanDiag and store diagonal clusters for sparse kmer!" << endl;

		flipCluster(sparse_clusts, fopts, liopt, 0);
		// if (sparse_clusts.size() == 0) {
		// uint32_t len = samples[a].genome->getLen(samples[a].idx);
		// sparse_clusts.push_back(cluster(0, len, 0, len, 0));
		// }

		// insert original clusters boundaries to bps
		for (auto&it : sparse_clusts) {
			bps.insert(it.yStart);
			bps.insert(it.yEnd);
			bps.insert(it.xStart);
			bps.insert(it.xEnd);
		}

		// insert the original breakpoints
		if (samples[a].left_bps.size() > 0) {
			for (auto&it : samples[a].left_bps)
				bps.insert(it);			
		}

		if (samples[a].right_bps.size() > 0) {
			for (auto&it : samples[a].right_bps)
				bps.insert(it);			
		}

		samples[a].hyperbreakpoints.clear();

		if (bps.size() > 0) {
		  	for (auto&it : bps) 
		    	samples[a].hyperbreakpoints.push_back(it);			
		}
	  	bps.clear();
	}
}

void sample::dump (string * readname, const fragopt_t &fopts)
{
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

void sample::substractClusters (bool dense)
{
	int size;

	if (dense) {
		size = dense_clusts.size();
		if (size > 0) 
			dense_clusts.resize( size / 2);
	}
	else {
		size = sparse_clusts.size();
		if (size > 0)
			sparse_clusts.resize(size / 2);
	}
}

// this function takes a bunch of breakpoints: [s, e) in bps
// returns cluster pivots: left_bps and right_bps
void Cluster_helper (uint32_t s, uint32_t e, const vector<uint32_t> &bps, vector<uint32_t> &trimInfo, vector<uint32_t> &right, vector<uint32_t> &left, const fragopt_t &fopts)
{
	uint32_t i;
	uint32_t prev_bp;

	assert(e >= s);
	assert(s < trimInfo.size() and e <= trimInfo.size());
	if (e == s + 1 or e == s) return;

	// find the rightmost cluster pivot
	prev_bp = bps[e - 1];
	i = e - 2; 
	while (i >= s) {
		if (bps[i] + fopts.clusterTrimedge >= prev_bp)
			trimInfo[i] = trimInfo[i + 1];
		else
			trimInfo[i] = i;
		prev_bp = bps[i];
		if (i == s) break;
		i--;
	}

	for (i = s; i < e; ++i) {
		if (trimInfo[i] == i)
			right.push_back(bps[i]);
	}

	// find the leftmost cluster pivot
	for (i = s; i < e; ++i) 
		trimInfo[i] = i;

	prev_bp = bps[s];
	i = s + 1; 
	while (i < e) {
		if (bps[i] <= prev_bp + fopts.clusterTrimedge)
			trimInfo[i] = trimInfo[i - 1];
		else
			trimInfo[i] = i;
		prev_bp = bps[i];
		i++;
	}

	for (i = s; i < e; ++i) {
		if (trimInfo[i] == i)
			left.push_back(bps[i]);
	}

	assert(right.size() == left.size());
}

// cluster breakpoints in 2 ways
// keep the rightmost cluster pivot in right_bps
// keep the leftmost cluster pivot in left_bps
void sample::clusterBreakpoints (const fragopt_t &fopts)
{
	// vector<uint32_t> clusterPivots;
	uint32_t sz = breakpoints.size();
	assert(sz > 0);
	vector<uint32_t> trimInfo(sz, 0);
	iota(trimInfo.begin(), trimInfo.end(), 0);
	// int i;
	uint32_t prev;

	assert(sz > 0);
	right_bps.clear(); left_bps.clear();
	Cluster_helper (0, sz, breakpoints, trimInfo, right_bps, left_bps, fopts);
	assert(right_bps.size() == left_bps.size());
	// breakpoints.clear();

/*
	right_bps.clear();
	// find the rightmost cluster pivot
	prev = breakpoints.back();
	i = breakpoints.size() - 2; 
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
			right_bps.push_back(breakpoints[i]);
	}

	// right_bps.clear();
	// right_bps = clusterPivots;

	// find the leftmost cluster pivot
	// clusterPivots.clear();
	left_bps.clear();
	iota(trimInfo.begin(), trimInfo.end(), 0);
	prev = breakpoints[0];
	i = 1; 
	while (i < breakpoints.size()) {
		if (breakpoints[i] <= prev + fopts.clusterTrimedge)
			trimInfo[i] = trimInfo[i - 1];
		else
			trimInfo[i] = i;
		prev = breakpoints[i];
		i++;
	}

	for (i = 0; i < trimInfo.size(); ++i) {
		if (trimInfo[i] == i)
			left_bps.push_back(breakpoints[i]);
	}

	// left_bps.clear();
	// left_bps = clusterPivots;
	// clusterPivots.clear();	

	assert(right_bps.size() == left_bps.size());
	// breakpoints.clear();
	*/
}

void sample::clusterBreakpointswithinInitalBreakpoints (const fragopt_t &fopts)
{
	// vector<uint32_t> clusterPivots;
	assert(hyperbreakpoints.size() > 0);
	vector<bool> remove(hyperbreakpoints.size(), 0);
	vector<uint32_t> right_idx(right_bps.size(), 0);
	vector<uint32_t> left_idx(left_bps.size(), 0);

 	right_hyperbps.clear();
 	left_hyperbps.clear();

	int i, j, k;
	j = 0;
	uint32_t prev;

	//delete every breakpoints between left_bps[i] and right_bps[i]
	//-----left0---right0--------------------------left1---right1----------
	j = 0;
	for (i = 0, j = 0; i < right_bps.size() and j < hyperbreakpoints.size(); ++i) {
		// [left_bps[i], right_bps[i]]
		while (hyperbreakpoints[j] <= left_bps[i]) { j++; }
		assert(hyperbreakpoints[j] > left_bps[i]);

		while (hyperbreakpoints[j] > left_bps[i] and hyperbreakpoints[j] < right_bps[i]) {
			remove[j] = 1;
			j++;
		}
		assert(hyperbreakpoints[j] == right_bps[i]);
	}

	uint32_t c = 0;
	for (i = 0; i < remove.size(); ++i) {
		if (remove[i] == 0) {
			hyperbreakpoints[c] = hyperbreakpoints[i];
			c++;
		}
	}
	hyperbreakpoints.resize(c);
	assert(hyperbreakpoints.size() > 0);


	vector<uint32_t> trimInfo(hyperbreakpoints.size());
	iota(trimInfo.begin(), trimInfo.end(), 0);

	// expand left_bps
	j = hyperbreakpoints.size() - 1;
	for (i = left_bps.size() - 1; i >= 0; --i) {
		while (hyperbreakpoints[j] > left_bps[i]) { j--; }
		assert(hyperbreakpoints[j] == left_bps[i]);
		left_idx[i] = j;

		prev = left_bps[i];
		k = j - 1;
		while (k >= 0 and hyperbreakpoints[k] > right_bps[i - 1]) {
			if (hyperbreakpoints[k] + fopts.clusterTrimedge >= prev) {
				trimInfo[k] = trimInfo[k + 1];
				left_bps[i] = hyperbreakpoints[k]; // update the left_bps[i]
				left_idx[i] = k;
			}
			else 
				break;
			prev = hyperbreakpoints[k];
			k--;
		}
		j = k;
	}

	// expand right_bps
	j = 0;
	for (i = 0; i < right_bps.size(); ++i) {
		while (hyperbreakpoints[j] < right_bps[i]) { j++; }
		assert(hyperbreakpoints[j] == right_bps[i]);
		right_idx[i] = j;

		prev = right_bps[i];
		k = j + 1;
		while (k < hyperbreakpoints.size() and ( i == right_bps.size() - 1 or (i < right_bps.size() - 1 and hyperbreakpoints[k] < left_bps[i + 1]) ) ) {
			if (hyperbreakpoints[k] <= prev + fopts.clusterTrimedge) {
				trimInfo[k] = trimInfo[k - 1];
				right_bps[i] = hyperbreakpoints[k];
				right_idx[i] = k;
			}
			else 
				break;
			prev = hyperbreakpoints[k];
			k++;
		}
		j = k;
	}

	// cluster the rest breakpoints
	// -----left0---right0--------------------------left1---right1----------
	for (i = 0; i < left_idx.size(); ++i) {
		if (i == 0) { 
			if (left_idx[0] > 0)
				Cluster_helper(0, left_idx[0], hyperbreakpoints, trimInfo, right_hyperbps, left_hyperbps, fopts);
		}
		else {
			Cluster_helper(right_idx[i - 1] + 1, left_idx[i], hyperbreakpoints, trimInfo, right_hyperbps, left_hyperbps, fopts);
			if (i == right_idx.size() - 1) {
				if (hyperbreakpoints.size() > right_idx.back() + 1)
					Cluster_helper(right_idx.back() + 1, hyperbreakpoints.size(), hyperbreakpoints, trimInfo, right_hyperbps, left_hyperbps, fopts);
			}
		}

		if (i < left_idx.size()) {
			left_hyperbps.push_back(left_bps[i]);
			right_hyperbps.push_back(right_bps[i]);
		}
	}


	// // cluster to the right
	// end = 0; prev_end = 0;
	// for (i = 0, j = 0; i < right_bps.size() and j < hyperbreakpoints.size(); ++i) {

	// 	if (j == 0) {
	// 		while (hyperbreakpoints[j] < right_bps[i]) { j++; }
	// 		assert(hyperbreakpoints[j] == right_bps[i]);
	// 		end = j;			
	// 	}

	// 	// while (hyperbreakpoints[j] < right_bps[i]) { j++; }
	// 	// assert(hyperbreakpoints[j] == right_bps[i]);
	// 	// right = j;

	// 	// (prev_end, end]
	// 	prev = hyperbreakpoints[end];
	// 	k = end - 1;
	// 	while (k > prev_end) {
	// 		if (hyperbreakpoints[k] + fopts.clusterTrimedge >= prev)
	// 			trimInfo[k] = trimInfo[k + 1];
	// 		else
	// 			trimInfo[k] = k;
	// 		prev = hyperbreakpoints[k];
	// 		k--;			
	// 	}
	// 	prev_end = end;
	// } 

	// for (i = 0; i < trimInfo.size(); ++i) {
	// 	if (trimInfo[i] == i)
	// 		clusterPivots.push_back(hyperbreakpoints[i]);
	// }

	// right_hyperbps.clear();
	// right_hyperbps = clusterPivots; // update breakpoints to cleaner version
	// // hyperbreakpoints = clusterPivots; 

	// // cluster to the left
	// trimInfo.clear();
	// clusterPivots.clear();
	// trimInfo.resize(hyperbreakpoints.size());
	// iota(trimInfo.begin(), trimInfo.end(), 0);
	// for (i = 1, j = 0; i <= left_bps.size() and j < hyperbreakpoints.size(); ++i) {

	// 	if (j == 0) {
	// 		while (hyperbreakpoints[j] < left_bps[i - 1]) { j++; }
	// 		assert(hyperbreakpoints[j] == left_bps[i - 1]);
	// 		left = j;			
	// 	}

	// 	if (i < left_bps.size()) {
	// 		while (hyperbreakpoints[j] < left_bps[i]) { j++; }
	// 		assert(hyperbreakpoints[j] == left_bps[i]);
	// 		right = j;			
	// 	}
	// 	else 
	// 		right = hyperbreakpoints.size();
		
	// 	// [left, right)
	// 	prev = hyperbreakpoints[left];
	// 	k = left + 1;
	// 	while (k < right) {
	// 		if (hyperbreakpoints[k] <= prev + fopts.clusterTrimedge)
	// 			trimInfo[k] = trimInfo[k - 1];
	// 		else
	// 			trimInfo[k] = k;
	// 		prev = hyperbreakpoints[k];
	// 		k++;			
	// 	}

	// 	left = right;
	// } 

	// for (i = 0; i < trimInfo.size(); ++i) {
	// 	if (trimInfo[i] == i)
	// 		clusterPivots.push_back(hyperbreakpoints[i]);
	// }

	// left_hyperbps.clear();
	// left_hyperbps = clusterPivots; // update breakpoints to cleaner version
	
	hyperbreakpoints.clear();
	assert(left_hyperbps.size() == right_hyperbps.size());
}

// void sample::clusterBreakpointswithinInitalBreakpoints (const fragopt_t &fopts, vector<uint32_t> &clusterPivots)
// {
// 	vector<uint32_t> trimInfo(hyperbreakpoints.size());
// 	iota(trimInfo.begin(), trimInfo.end(), 0);

// 	int i, j, k;
// 	int left, right;
// 	j = 0;
// 	uint32_t prev;

// 	// two pointers
// 	// cluster to the right
// 	for (i = 1; i < breakpoints.size() and j < hyperbreakpoints.size(); ++i) {

// 		if (j == 0) {
// 			while (hyperbreakpoints[j] < breakpoints[i - 1]) { j++; }
// 			assert(hyperbreakpoints[j] == breakpoints[i - 1]);
// 			left = j;			
// 		}

// 		while (hyperbreakpoints[j] < breakpoints[i]) { j++; }
// 		assert(hyperbreakpoints[j] == breakpoints[i]);
// 		right = j;

// 		// (left, right]
// 		prev = hyperbreakpoints[right];
// 		k = right - 1;
// 		while (k > left) {
// 			if (hyperbreakpoints[k] + fopts.clusterTrimedge >= prev)
// 				trimInfo[k] = trimInfo[k + 1];
// 			else
// 				trimInfo[k] = k;
// 			prev = hyperbreakpoints[k];
// 			k--;			
// 		}

// 		left = right;
// 	} 

// 	for (i = 0; i < trimInfo.size(); ++i) {
// 		if (trimInfo[i] == i)
// 			clusterPivots.push_back(hyperbreakpoints[i]);
// 	}

// 	hyperbreakpoints = clusterPivots; // update breakpoints to cleaner version


// 	// cluster to the left
// 	trimInfo.clear();
// 	clusterPivots.clear();
// 	trimInfo.resize(hyperbreakpoints.size());
// 	iota(trimInfo.begin(), trimInfo.end(), 0);
// 	for (i = 1, j = 0; i <= breakpoints.size() and j < hyperbreakpoints.size(); ++i) {

// 		if (j == 0) {
// 			while (hyperbreakpoints[j] < breakpoints[i - 1]) { j++; }
// 			assert(hyperbreakpoints[j] == breakpoints[i - 1]);
// 			left = j;			
// 		}

// 		if (i < breakpoints.size()) {
// 			while (hyperbreakpoints[j] < breakpoints[i]) { j++; }
// 			assert(hyperbreakpoints[j] == breakpoints[i]);
// 			right = j;			
// 		}
// 		else 
// 			right = hyperbreakpoints.size();
		
// 		// [left, right)
// 		prev = hyperbreakpoints[left];
// 		k = left + 1;
// 		while (k < right) {
// 			if (hyperbreakpoints[k] <= prev + fopts.clusterTrimedge)
// 				trimInfo[k] = trimInfo[k - 1];
// 			else
// 				trimInfo[k] = k;
// 			prev = hyperbreakpoints[k];
// 			k++;			
// 		}

// 		left = right;
// 	} 

// 	for (i = 0; i < trimInfo.size(); ++i) {
// 		if (trimInfo[i] == i)
// 			clusterPivots.push_back(hyperbreakpoints[i]);
// 	}

// 	hyperbreakpoints = clusterPivots; // update breakpoints to cleaner version
// }

bool Maxoverlap (uint32_t b_s, uint32_t b_e, uint32_t a_s, uint32_t a_e, double frac)
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

bool Minoverlap (uint32_t b_s, uint32_t b_e, uint32_t a_s, uint32_t a_e, double frac)
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
		if ( min(ovp / denomA, ovp / denomB) >= frac) 
			return true;
		else 
			return false;
}
void addEdge(const fragopt_t &fopts, const vector<cluster> &clusts, vector<uint32_t> &left_bps, vector<uint32_t> &right_bps, vector<vector<uint32_t>> &adlist, uint32_t s)
{
	uint32_t i;
	vector<uint32_t> xhits, yhits;

	for (auto&it : clusts) {
		xhits.clear();
		yhits.clear();
		for (i = 1; i < left_bps.size(); ++i) { // larger interval: [left_bps[i - 1], right_bps[i]]; smaller interval: [right_bps[i-1], left_bps[i]]
			// if (Maxoverlap(right_bps[i - 1], left_bps[i], it.xStart, it.xEnd, fopts.ovpfrac))
			// 	xhits.push_back(i - 1);
			// if (Maxoverlap(right_bps[i - 1], left_bps[i], it.yStart, it.yEnd, fopts.ovpfrac))
			// 	yhits.push_back(i - 1);
			if (Maxoverlap(left_bps[i - 1], right_bps[i], it.xStart, it.xEnd, fopts.ovpfrac) or Maxoverlap(right_bps[i - 1], left_bps[i], it.xStart, it.xEnd, fopts.ovpfrac)
				or Maxoverlap(left_bps[i - 1], left_bps[i], it.xStart, it.xEnd, fopts.ovpfrac) or Maxoverlap(right_bps[i - 1], right_bps[i], it.xStart, it.xEnd, fopts.ovpfrac))
				xhits.push_back(i - 1);
			if (Maxoverlap(left_bps[i - 1], right_bps[i], it.yStart, it.yEnd, fopts.ovpfrac) or Maxoverlap(right_bps[i - 1], left_bps[i], it.yStart, it.yEnd, fopts.ovpfrac)
				or Maxoverlap(left_bps[i - 1], left_bps[i], it.yStart, it.yEnd, fopts.ovpfrac) or Maxoverlap(right_bps[i - 1], right_bps[i], it.yStart, it.yEnd, fopts.ovpfrac))
				yhits.push_back(i - 1);
		}
		assert(xhits.size() == yhits.size());

		for (i = 0; i < xhits.size(); ++i) {
			adlist[xhits[i] + s].push_back(yhits[i] + s);
			adlist[yhits[i] + s].push_back(xhits[i] + s);
		}
	}
}

void sample::unifyIntv(const fragopt_t &fopts, graph &superGraph, bool dense)
{
	if (dense)
		addEdge(fopts, dense_clusts, left_bps, right_bps, superGraph.adlist_dense, s);
	else
		addEdge(fopts, sparse_clusts, left_hyperbps, right_hyperbps, superGraph.adlist_sparse, s);
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

void acrosample::across_process (sample &sample_a, sample &sample_b, const uint32_t a, const uint32_t b, const idx_t &mi_s, const idx_t &mi_l, 
		const fragopt_t &fopts, const idxopt_t &siopt, const idxopt_t &liopt, bool self, bool dense)
{
	cerr << "process across sample " << a << " " << b << endl;
	if (dense) {
		// get diagonal cluster matches
		vector<hit> dense_fhits, dense_rhits;
		rf_hit(sample_a, sample_b, a, b, mi_s, dense_fhits, dense_rhits, fopts, siopt, 0); // find rhits also
		cerr << "finish rf_hit for dense kmer!" << endl;

		// clusters dense_clusts;
		cleanDiag(dense_fhits, dense_clusts, fopts, siopt, 0);
		cleanDiag(dense_rhits, dense_clusts, fopts, siopt, 1);
		assert(dense_clusts.size() > 0);

		// if (fopts.debug) {
			ofstream clust("pair_hits.bed", ios_base::app);
			for (int m = 0; m < dense_fhits.size(); m++) {
				clust << dense_fhits[m].x << "\t" << dense_fhits[m].y << "\t" << dense_fhits[m].x + mi_s.k << "\t"
					   << dense_fhits[m].y + mi_s.k << "\t" << mi_s.k << "\t"  
					   << "0"  << "\t" << "1" << "\t" << *(sample_a.readname) + " " + *(sample_b.readname) << endl;				
			}
			for (int m = 0; m < dense_rhits.size(); m++) {
				clust << dense_rhits[m].x << "\t" << dense_rhits[m].y << "\t" << dense_rhits[m].x + mi_s.k << "\t"
					   << dense_rhits[m].y + mi_s.k << "\t" << mi_s.k << "\t"  
					   << "1" << "\t" << "1" << "\t" << *(sample_a.readname) + " " + *(sample_b.readname) << endl;				
			}
			clust.close();			
		// }

		dense_fhits.clear();
		dense_rhits.clear();
		cerr << "finish cleanDiag and store diagonal clusters for dense kmer!" << endl;		
	}
	else {
		vector<hit> sparse_fhits, sparse_rhits;
		rf_hit(sample_a, sample_b, a, b, mi_l, sparse_fhits, sparse_rhits, fopts, liopt, 0); // find rhits also
		cerr << "finish rf_hit for sparse kmer!" << endl;

		// clusters sparse_clusts;
		cleanDiag(sparse_fhits, sparse_clusts, fopts, liopt, 0);
		cleanDiag(sparse_rhits, sparse_clusts, fopts, liopt, 1);
		assert(sparse_clusts.size() > 0);

		// if (fopts.debug) {
			ofstream clust("pair_hits.bed", ios_base::app);
			for (int m = 0; m < sparse_fhits.size(); m++) {
				clust << sparse_fhits[m].x << "\t" << sparse_fhits[m].y << "\t" << sparse_fhits[m].x + mi_s.k << "\t"
					   << sparse_fhits[m].y + mi_l.k << "\t" << mi_l.k << "\t"  
					   << "0"  << "\t" << "0" << "\t" << *(sample_a.readname) + " " + *(sample_b.readname) << endl;				
			}
			for (int m = 0; m < sparse_rhits.size(); m++) {
				clust << sparse_rhits[m].x << "\t" << sparse_rhits[m].y << "\t" << sparse_rhits[m].x + mi_s.k << "\t"
					   << sparse_rhits[m].y + mi_l.k << "\t" << mi_l.k << "\t"  
					   << "1" << "\t" << "0" << "\t" << *(sample_a.readname) + " " + *(sample_b.readname) << endl;				
			}
			clust.close();			
		// }

		sparse_fhits.clear();
		sparse_rhits.clear();
		cerr << "finish cleanDiag and store diagonal clusters for sparse kmer!" << endl;		
	}
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

void trimOneCluster (cluster &clust, uint32_t left, uint32_t right, bool axis = 1) 
{
	if (axis == 1) { // x-axis
		if (clust.xStart < left)
			clust.xStart = left;
		if (clust.xEnd > right)
			clust.xEnd = right;
	}
	else {
		if (clust.yStart < left)
			clust.yStart = left;
		if (clust.yEnd > right)
			clust.yEnd = right;
	}
}

void acrosample::trimclusters(bool dense)
{
	uint32_t left, right;
	if (dense) {
		left = pi->breakpoints[0];
		right = pi->breakpoints.back();
		for (auto&it : dense_clusts) {
			trimOneCluster(it, left, right);
		}
		left = pj->breakpoints[0];
		right = pj->breakpoints.back();		
		for (auto&it : dense_clusts) {
			trimOneCluster(it, left, right, 0);
		}
	}
	else {
		left = pi->breakpoints[0];
		right = pi->breakpoints.back();
		for (auto&it : sparse_clusts) {
			trimOneCluster(it, left, right);
		}
		left = pj->breakpoints[0];
		right = pj->breakpoints.back();		
		for (auto&it : sparse_clusts) {
			trimOneCluster(it, left, right, 0);
		}		
	}
}

void acrosample::dump (string * readname_i, string * readname_j, const fragopt_t &fopts, uint32_t len, bool dense)
{
	ofstream pclust("pair_cluster.bed", ios_base::app);
	for (int m = 0; m < dense_clusts.size(); m++) {
		pclust << dense_clusts[m].xStart << "\t" << dense_clusts[m].yStart << "\t" << dense_clusts[m].xEnd << "\t"
			   << dense_clusts[m].yEnd << "\t" << dense_clusts[m].xEnd - dense_clusts[m].xStart << "\t"  
			   << dense_clusts[m].strand  << "\t" << m << "\t" << "1" << "\t" << *readname_i + " " + *readname_j << endl;				
	}
	for (int m = 0; m < sparse_clusts.size(); m++) {
		pclust << sparse_clusts[m].xStart << "\t" << sparse_clusts[m].yStart << "\t" << sparse_clusts[m].xEnd << "\t"
			   << sparse_clusts[m].yEnd << "\t" << sparse_clusts[m].xEnd - sparse_clusts[m].xStart << "\t"  
			   << sparse_clusts[m].strand  << "\t" << m << "\t" << "0" << "\t" << *readname_i + " " + *readname_j << endl;				
	}
	pclust.close();				

	ofstream clust("cluster_acrosssample.bed", ios_base::app);
	// cerr << "across sample dense_clusts.size(): " << dense_clusts.size() << endl;
	// cerr << "across sample sparse_clusts.size(): " << sparse_clusts.size() << endl;

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
	if (dense) {
	  	for (auto& it : pi->breakpoints)
	    	fclust << it << "\t" << *readname_i << "\t" << *readname_j << endl;	
	  	for (auto& it : pj->breakpoints)
	    	fclust << it << "\t" << *readname_j << "\t" << *readname_i << endl;
		fclust.close();				
	}
	else {
	  	for (auto& it : pi->hyperbreakpoints)
	    	fclust << it << "\t" << *readname_i << "\t" << *readname_j << endl;	
	  	for (auto& it : pj->hyperbreakpoints)
	    	fclust << it << "\t" << *readname_j << "\t" << *readname_i << endl;
		fclust.close();				
	}
}

void addEdge_acros (const fragopt_t &fopts, const vector<cluster> &clusts, sample *pi, sample *pj, vector<vector<uint32_t>> &adlist, 
					const vector<uint32_t> &pi_leftbps, const vector<uint32_t> &pi_rightbps, const vector<uint32_t> &pj_leftbps, const vector<uint32_t> &pj_rightbps)
{
	uint32_t c;
	uint32_t t;
	vector<uint32_t> xhits, yhits;
	uint32_t m, n;

	for (auto&it : clusts) {
		xhits.clear();
		yhits.clear();
		// for (c = 1; c < pi_leftbps.size(); ++c) {
		// 	if (Maxoverlap(pi_rightbps[c - 1], pi_leftbps[c], it.xStart, it.xEnd, fopts.ovpfrac))
		// 		xhits.push_back(c - 1);
		// }

		// for (c = 1; c < pj_leftbps.size(); ++c) {
		// 	if (Maxoverlap(pj_rightbps[c - 1], pj_leftbps[c], it.yStart, it.yEnd, fopts.ovpfrac))
		// 		yhits.push_back(c - 1);
		// }
		for (c = 1; c < pi_leftbps.size(); ++c) {
			if (Maxoverlap(pi_leftbps[c - 1], pi_rightbps[c], it.xStart, it.xEnd, fopts.ovpfrac) or Maxoverlap(pi_rightbps[c - 1], pi_leftbps[c], it.xStart, it.xEnd, fopts.ovpfrac)
				or Maxoverlap(pi_leftbps[c - 1], pi_leftbps[c], it.xStart, it.xEnd, fopts.ovpfrac) or Maxoverlap(pi_rightbps[c - 1], pi_rightbps[c], it.xStart, it.xEnd, fopts.ovpfrac))
				xhits.push_back(c - 1);
		}

		for (c = 1; c < pj_leftbps.size(); ++c) {
			if (Maxoverlap(pj_leftbps[c - 1], pj_rightbps[c], it.yStart, it.yEnd, fopts.ovpfrac) or Maxoverlap(pj_rightbps[c - 1], pj_leftbps[c], it.yStart, it.yEnd, fopts.ovpfrac)
				or Maxoverlap(pj_leftbps[c - 1], pj_leftbps[c], it.yStart, it.yEnd, fopts.ovpfrac) or Maxoverlap(pj_rightbps[c - 1], pj_rightbps[c], it.yStart, it.yEnd, fopts.ovpfrac))
				yhits.push_back(c - 1);
		}

		assert(xhits.size() == yhits.size());

		if (it.strand == 0) {
			for (c = 0; c < xhits.size(); ++c) {
				m = xhits[c] + pi->s;
				n = yhits[c] + pj->s;
				adlist[m].push_back(n);
				adlist[n].push_back(m);
				// uf->merge(m + 2, n + 2); // same strand
			}
		}
		else {
			for (c = 0; c < xhits.size(); ++c) {
				t = xhits.size() - 1 - c;
				m = xhits[c] + pi->s;
				n = yhits[t] + pj->s;
				adlist[m].push_back(n);
				adlist[n].push_back(m);
			}			
		}
	}	
}

void acrosample::unifyIntv(const fragopt_t &fopts, graph &superGraph, bool dense)
{
	if (dense)
		addEdge_acros(fopts, dense_clusts, pi, pj, superGraph.adlist_dense, pi->left_bps, pi->right_bps, pj->left_bps, pj->right_bps);
	else
		addEdge_acros(fopts, sparse_clusts, pi, pj, superGraph.adlist_sparse, pi->left_hyperbps, pi->right_hyperbps, pj->left_hyperbps, pj->right_hyperbps);
}

void dumpGraph (vector<sample> &samples, graph &G, const fragopt_t &fopts, bool dense)
{
	// if (!fopts.debug) return;
	uint32_t i, k;

	if (dense) {
		ofstream fclust("dense_assignment.bed", ios_base::app);
		for (i = 0; i < samples.size(); ++i) {
			fclust << *(samples[i].readname) << ": "; 
			for (k = samples[i].s; k < samples[i].e; ++k) {
				if (k != samples[i].e - 1)
					fclust << G.colors_dense[k] << ", ";
				else 
					fclust << G.colors_dense[k] << endl;			
			}
		}
		fclust.close();		

		ofstream lclust("dense_lens.bed", ios_base::app);
		for (i = 0; i < samples.size(); ++i) {
			lclust << *(samples[i].readname) << ": "; 
			for (k = samples[i].s; k < samples[i].e; ++k) {
				if (k != samples[i].e - 1) {
					lclust << G.lens[k] << ", ";

				}
				else {
					lclust << G.lens[k] << endl;
				}
			}
		}
		lclust.close();			
	}
	else {
		ofstream hclust("sparse_assignment.bed", ios_base::app);
		for (i = 0; i < samples.size(); ++i) {
			hclust << *(samples[i].readname) << ": "; 
			for (k = samples[i].s; k < samples[i].e; ++k) {
				if (k != samples[i].e - 1)
					hclust << G.colors_sparse[k] << ", ";
				else 
					hclust << G.colors_sparse[k] << endl;			
			}
		}
		hclust.close();	

		ofstream slclust("sparse_lens.bed", ios_base::app);
		for (i = 0; i < samples.size(); ++i) {
			slclust << *(samples[i].readname) << ": "; 
			for (k = samples[i].s; k < samples[i].e; ++k) {
				if (k != samples[i].e - 1) {
					slclust << G.lens[k] << ", ";

				}
				else {
					slclust << G.lens[k] << endl;
				}
			}
		}
		slclust.close();		
	}

	// ofstream kclust("assignment.bed", ios_base::app);
	// for (i = 0; i < samples.size(); ++i) {
	// 	kclust << *(samples[i].readname) << ": "; 
	// 	for (k = samples[i].s; k < samples[i].e; ++k) {
	// 		if (k != samples[i].e - 1) {
	// 			if (G.colors_sparse[k] != 0)
	// 				kclust << G.colors_sparse[k] << ", ";
	// 			else
	// 				kclust << G.colors_dense[k] << ", ";
	// 		}
	// 		else {
	// 			if (G.colors_sparse[k] != 0)
	// 				kclust << G.colors_sparse[k] << endl;
	// 			else
	// 				kclust << G.colors_dense[k] << endl;
	// 		}
	// 	}
	// }
	// kclust.close();	

			
}


void FindBasicRepeatUnit(vector<sample> &samples, vector<vector<acrosample>> &acrosamples, Genome &genome, const idx_t &mi_s, const idx_t &mi_l, 
						const idxopt_t &siopt, const idxopt_t &liopt, const fragopt_t &fopts)
{
	uint32_t i, j;
	uint32_t n = genome.getSize();
	vector<uint32_t> clusterPivots;

	// Get the breakpoints for each sample dense_clusts
	for (i = 0; i < n; ++i) {
		cerr << "process sample: " << *genome.getName(i) << " " << i << endl;
		samples[i].init(i, genome.getName(i), &genome);
		samples[i].process(samples, i, i, mi_s, mi_l, fopts, siopt, liopt, 1, 1); // dense
		// samples[i].dump(genome.getName(i), fopts);
	}

	// get the dense_clusts for across-sample 
	for (i = 0; i < n; ++i) {
		acrosamples[i].resize(n - i);
		for (j = i + 1; j < n; ++j) {
			// assert(samples[i].breakpoints.size() > 0 and samples[j].breakpoints.size() > 0);
			acrosamples[i][j - i - 1].init(i, j, &samples[i], &samples[j]);
			acrosamples[i][j - i - 1].across_process(samples[i], samples[j], i, j, mi_s, mi_l, fopts, siopt, liopt, 0, 1); // dense
			acrosamples[i][j - i - 1].decideStrand();
			acrosamples[i][j - i - 1].dump(genome.getName(i), genome.getName(j), fopts, genome.getLen(j), 1);
		}
	}	

	// project every sample to samples[i]
	if (n > 1) {
		for (i = 0; i < n; ++i) {
			for (j = 0; j < n; ++j) {
				// project samples[j] to samples[i]
				if (i == j) continue;
				if (j > i) { // samples[j] is y-axis
					project(samples[j], samples[i], acrosamples[i][j - i - 1].dense_clusts, acrosamples[i][j - i - 1].sparse_clusts, fopts, 0, 1);
				}
				else { // samples[j] is x-axis
					project(samples[j], samples[i], acrosamples[j][i - j - 1].dense_clusts, acrosamples[j][i - j - 1].sparse_clusts, fopts, 1, 1);
				}
			}
		}		
	}

	// for (i = 0; i < n; ++i) {
	// 	samples[i].clusterBreakpoints(fopts);
	// }
	
	// project samples[i] back to other samples
	if (n > 1) {
		for (i = 0; i < n; ++i) {
			for (j = 0; j < n; ++j) {
				// project samples[i] to samples[j]
				if (i == j) continue;
				if (j < i) { // samples[i] is y-axis
					project(samples[i], samples[j], acrosamples[j][i - j - 1].dense_clusts, acrosamples[j][i - j - 1].sparse_clusts, fopts, 0, 1);
				}
				else { // samples[i] is x-axis
					project(samples[i], samples[j], acrosamples[i][j - i - 1].dense_clusts, acrosamples[i][j - i - 1].sparse_clusts, fopts, 1, 1);
				}
			}
		}		
	}

	// for (i = 0; i < n; ++i) {
	// 	// samples[i].dump(genome.getName(i), fopts);
	// 	for (j = i + 1; j < n; ++j) {
	// 		acrosamples[i][j - i - 1].dump(genome.getName(i), genome.getName(j), fopts, genome.getLen(j), 1);
	// 	}
	// }

	// for (i = 0; i < n; ++i) {
	// 	samples[i].clusterBreakpoints(fopts);
	// }

	// project self again
	for (i = 0; i < n; ++i) {
		assert(samples[i].breakpoints.size() > 0);
		selfproject (samples[i], samples[i].dense_clusts, samples[i].sparse_clusts, fopts, 1);
		samples[i].clusterBreakpoints(fopts);
	}

	cerr << "finish projecting breakpoints" << endl;
}

void FindsmallerRepeatUnit (vector<sample> &samples, vector<vector<acrosample>> &acrosamples, Genome &genome, const idx_t &mi_s, const idx_t &mi_l, 
						const idxopt_t &siopt, const idxopt_t &liopt, const fragopt_t &fopts)
{
	uint32_t i, j;
	uint32_t n = genome.getSize();
	vector<uint32_t> clusterPivots;

	// Get the hyperbreakpoints for each sample sparse_clusts
	for (i = 0; i < n; ++i) {
		cerr << "process sample: " << *genome.getName(i) << " " << i << endl;
		samples[i].process(samples, i, i, mi_s, mi_l, fopts, siopt, liopt, 1, 0); // sparse
		// samples[i].dump(genome.getName(i), fopts);
	}


	// get the sparse_clusts for across-sample 
	for (i = 0; i < n; ++i) {
		for (j = i + 1; j < n; ++j) {
			acrosamples[i][j - i - 1].across_process(samples[i], samples[j], i, j, mi_s, mi_l, fopts, siopt, liopt, 0, 0); // sparse
			acrosamples[i][j - i - 1].dump(genome.getName(i), genome.getName(j), fopts, genome.getLen(j), 0);
		}
	}	

	// project every sample to samples[i]
	if (n > 1) {
		for (i = 0; i < n; ++i) {
			for (j = 0; j < n; ++j) {
				// project samples[j] to samples[i]
				if (i == j) continue;
				if (j > i) { // samples[j] is y-axis
					project(samples[j], samples[i], acrosamples[i][j - i - 1].dense_clusts, acrosamples[i][j - i - 1].sparse_clusts, fopts, 0, 0);
				}
				else { // samples[j] is x-axis
					project(samples[j], samples[i], acrosamples[j][i - j - 1].dense_clusts, acrosamples[j][i - j - 1].sparse_clusts, fopts, 1, 0);
				}
			}
		}		
	}

	// for (i = 0; i < n; ++i) {
	// 	clusterPivots.clear();
	// 	samples[i].clusterBreakpointswithinInitalBreakpoints(fopts, clusterPivots);
	// }

	// project samples[i] back to other samples
	if (n > 1) {
		for (i = 0; i < n; ++i) {
			for (j = 0; j < n; ++j) {
				// project samples[i] to samples[j]
				if (i == j) continue;
				if (j < i) { // samples[i] is y-axis
					project(samples[i], samples[j], acrosamples[j][i - j - 1].dense_clusts, acrosamples[j][i - j - 1].sparse_clusts, fopts, 0, 0);
				}
				else { // samples[i] is x-axis
					project(samples[i], samples[j], acrosamples[i][j - i - 1].dense_clusts, acrosamples[i][j - i - 1].sparse_clusts, fopts, 1, 0);
				}
			}
		}		
	}

	for (i = 0; i < n; ++i) {
		assert(samples[i].hyperbreakpoints.size() > 0);
		// samples[i].clusterBreakpointswithinInitalBreakpoints(fopts, clusterPivots);
	}

	// project self again
	for (i = 0; i < n; ++i) {
		selfproject (samples[i], samples[i].dense_clusts, samples[i].sparse_clusts, fopts, 0);
		for (j = i + 1; j < n; ++j) {
			acrosamples[i][j - i - 1].dump(genome.getName(i), genome.getName(j), fopts, genome.getLen(j), 0);
		}
		samples[i].clusterBreakpointswithinInitalBreakpoints(fopts);
	}

	// if (fopts.debug) {
	// 	for (i = 0; i < n; ++i) {
	// 		samples[i].dump(genome.getName(i), fopts);
	// 		for (j = i + 1; j < n; ++j) {
	// 			acrosamples[i][j - i - 1].dump(genome.getName(i), genome.getName(j), fopts, genome.getLen(j), 0);
	// 		}
	// 	}
	// }

	cerr << "finish projecting breakpoints" << endl;

}

void graphConstruction (graph &G, uint32_t n, vector<sample> &samples, vector<vector<acrosample>> &acrosamples, const fragopt_t &fopts, bool dense)
{
	uint32_t i, j;
	for (i = 0; i < n; ++i) {
		// // merge breakpoints if adjacent breakpoints are within 500bp
		// mergeCloseBreakpoints(samples[i], dense)

		// substract self-self dotplot clusters
		samples[i].substractClusters(dense);
		if (dense)
			G.insertInvt(samples[i].left_bps, samples[i].right_bps, samples[i].idx, samples[i].s, samples[i].e);
		else
			G.insertInvt(samples[i].left_hyperbps, samples[i].right_hyperbps, samples[i].idx, samples[i].s, samples[i].e);
	}

	G.init();
	for (i = 0; i < n; ++i) 
		samples[i].unifyIntv(fopts, G, dense);

	for (i = 0; i < n; ++i) {
		for (j = i + 1; j < n; ++j) 
			acrosamples[i][j - i - 1].unifyIntv(fopts, G, dense);
	}
	cerr << "finish constructing the graph" << endl;

	// find connected components in the superGraph
	G.findConnetedComponents(dense);
	dumpGraph(samples, G, fopts, dense);	

}









