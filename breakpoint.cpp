#include "breakpoint.h"
#include "cluster.h"
#include "unionfind.h"
#include <numeric>
#include <vector>
#include <fstream>
#include <set>


void trimBreakpoints(const vector<uint32_t> &bps, vector<uint32_t> &trimInfo)
{
	uint32_t prev = bps.back();
	int i = bps.size() - 2; 
	while (i >= 0) {
		if (bps[i] + 200 >= prev)
			trimInfo[i] = trimInfo[i + 1];
		else
			trimInfo[i] = i;
		prev = bps[i];
		i--;
	}
}

void modifyClusterBoundaries (vector<cluster> &clusts, uint32_t original, uint32_t after, bool axis) 
{
	uint32_t c;
	for (c = 0; c < clusts.size(); ++c) {
		assert(after > original);
		if (axis == 0) { // y-axis
			if (clusts[c].yStart == original) {
				clusts[c].yStart = after;
				clusts[c].xStart += (after - original);
			}
			else if (clusts[c].yEnd == original) {
				clusts[c].yEnd = after;
				clusts[c].xEnd += (after - original);
			}			
		}
		else {
			if (clusts[c].xStart == original) {
				clusts[c].xStart = after;
				clusts[c].yStart += (after - original);
			}
			else if (clusts[c].xEnd == original) {
				clusts[c].xEnd = after;
				clusts[c].yEnd += (after - original);
			}
		}

	}
}

void trimClusters(string * readname, vector<uint32_t> &bps, vector<uint32_t> &trimInfo, vector<cluster> &dense_clusts, vector<cluster> &sparse_clusts, bool axis = 0)
{
	uint32_t i;
	vector<bool> remove(bps.size(), 0);
	for (i = 0; i < trimInfo.size(); ++i) {
		if (trimInfo[i] == i) continue;
		remove[i] = 1;
		modifyClusterBoundaries(dense_clusts, bps[i], bps[trimInfo[i]], axis);
		modifyClusterBoundaries(sparse_clusts, bps[i], bps[trimInfo[i]], axis);
	}

	uint32_t c = 0;
	for (i = 0; i < remove.size(); ++i) {
		if (remove[i] == 0) {
			bps[c] = bps[i];
			c++;
		}
	}
	bps.resize(c);
}

void trimOnY(string * readname, vector<cluster> &dense_clusts, vector<cluster> &sparse_clusts, vector<uint32_t> &bps, const fragopt_t &fopts)
{
	set<uint32_t> breakpoints;
	uint32_t j;
	for (j = 0; j < dense_clusts.size(); ++j) {
		breakpoints.insert(dense_clusts[j].yStart);
		breakpoints.insert(dense_clusts[j].yEnd);
	}
	for (j = 0; j < sparse_clusts.size(); ++j) {
		breakpoints.insert(sparse_clusts[j].yStart);
		breakpoints.insert(sparse_clusts[j].yEnd);			
	}
  	for (auto it = breakpoints.begin(); it != breakpoints.end(); ++it)
    	bps.push_back(*it);
	// if (fopts.debug) {
	// 	ofstream fclust("lines.bed");
	//   	for (auto it = breakpoints.begin(); it != breakpoints.end(); ++it)
	//     	fclust << *it << endl;
	// 	fclust.close();			
	// }
    breakpoints.clear();

	// trim the breakpoints on y-axis
	vector<uint32_t> trimInfo(bps.size());
	iota(trimInfo.begin(), trimInfo.end(), 0);
	trimBreakpoints(bps, trimInfo);
	trimClusters(readname, bps, trimInfo, dense_clusts, sparse_clusts, 0);
	if (fopts.debug) {
		ofstream fclust("trimlines_Y.bed", ios_base::app);
	  	for (auto& it : bps)
	    	fclust << it << "\t" << *readname << endl;
		fclust.close();			
	}	
}

void trimOnX(string * readname, vector<cluster> &dense_clusts, vector<cluster> &sparse_clusts, vector<uint32_t> &bps_Y, vector<uint32_t> &bps_X, const fragopt_t &fopts)
{
	set<uint32_t> breakpoints_Y;
	set<uint32_t> breakpoints_X;
	uint32_t i, diff;
  	for (i = 0; i != bps_Y.size(); ++i)
    	breakpoints_Y.insert(bps_Y[i]);

	for (i = 0; i < dense_clusts.size(); ++i) {
		auto its = breakpoints_Y.lower_bound(dense_clusts[i].yStart); 
		auto ite = breakpoints_Y.lower_bound(dense_clusts[i].yEnd);
		assert(*its == dense_clusts[i].yStart and *ite == dense_clusts[i].yEnd);

		for (auto it = its; it != ite; ++it)
			diff = (*it > *its) ? (*it - *its) : 0;
			breakpoints_X.insert(diff + dense_clusts[i].xStart);
	}
	for (i = 0; i < sparse_clusts.size(); ++i) {
		auto its = breakpoints_Y.lower_bound(sparse_clusts[i].yStart); 
		auto ite = breakpoints_Y.lower_bound(sparse_clusts[i].yEnd);
		assert(*its == sparse_clusts[i].yStart and *ite == sparse_clusts[i].yEnd);

		for (auto it = its; it != ite; ++it)
			diff = (*it > *its) ? (*it - *its) : 0;
			breakpoints_X.insert(diff + sparse_clusts[i].xStart);	
	}
  	for (i = 0; i != bps_Y.size(); ++i)
    	breakpoints_X.insert(bps_Y[i]);

  	for (auto it = breakpoints_X.begin(); it != breakpoints_X.end(); ++it)
    	bps_X.push_back(*it);

    breakpoints_X.clear();
    breakpoints_Y.clear();

	// trim the breakpoints on x-axis
	vector<uint32_t> trimInfo(bps_X.size());
	iota(trimInfo.begin(), trimInfo.end(), 0);
	trimBreakpoints(bps_X, trimInfo);
	trimClusters(readname, bps_X, trimInfo, dense_clusts, sparse_clusts, 1);
	cerr << " get all the breakpoints on y, x-axis!" << endl;

	// if (fopts.debug) {
	// 	ofstream fclust("trimlines_X.bed");
	//   	for (auto& it : bps_X)
	//     	fclust << it << endl;
	// 	fclust.close();			
	// }	

}










