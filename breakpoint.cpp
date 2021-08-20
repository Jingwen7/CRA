#include "breakpoint.h"
#include "rfpriv.h"
#include "cluster.h"
#include "unionfind.h"
#include <numeric>
#include <vector>
#include <fstream>
#include <set>

void checkBreakpoints_Clusters (vector<uint32_t> &bps, vector<cluster> &clusts, bool axis) 
{
	uint32_t i;
	set<uint32_t> s;
	for (i = 0; i < bps.size(); ++i) {
		s.insert(bps[i]);
	}
	for (i = 0; i < clusts.size(); ++i) {
		if (axis == 0) {
			assert (s.count(clusts[i].yStart) != 0);
			assert(s.count(clusts[i].yEnd) != 0);
		}
		else {
			assert (s.count(clusts[i].xStart) != 0);
			assert(s.count(clusts[i].xEnd) != 0);	
		}
	}
}

void checkBreakpoints_Clusters (set<uint32_t> &s, vector<cluster> &clusts, bool axis) 
{
	uint32_t i;
	for (i = 0; i < clusts.size(); ++i) {
		if (axis == 0) {
			assert (s.count(clusts[i].yStart) != 0);
			assert(s.count(clusts[i].yEnd) != 0);
		}
		else {
			assert (s.count(clusts[i].xStart) != 0);
			assert(s.count(clusts[i].xEnd) != 0);	
		}
	}
}

void insertPoint (vector<cluster> &clusts, set<uint32_t> &s, bool axis) 
{
	uint32_t i;
	for (i = 0; i < clusts.size(); ++i) {
		if (axis == 0) {
			s.insert(clusts[i].xStart);
			s.insert(clusts[i].xEnd);
		}
		else{
			s.insert(clusts[i].yStart);
			s.insert(clusts[i].yEnd);
		}
	}
}

void insertPoint (vector<cluster> &clusts, vector<uint32_t> &s, bool axis) 
{
	uint32_t i;
	for (i = 0; i < clusts.size(); ++i) {
		if (axis == 0) {
			s.push_back(clusts[i].xStart);
			s.push_back(clusts[i].xEnd);
		}
		else{
			s.push_back(clusts[i].yStart);
			s.push_back(clusts[i].yEnd);
		}
	}
}

template<typename T>
void REsize(vector<T> &original, const vector<bool> &remove)
{
	uint32_t c = 0; 
	uint32_t i;
	for (i = 0; i < remove.size(); ++i) {
		if (remove[i] == 0) {
			original[c] = original[i];
			c++;
		}
	}
	original.resize(c);	
}

void clusterBreakpoints (const vector<uint32_t> &bps, vector<uint32_t> &trimInfo, const fragopt_t &fopts)
{
	uint32_t prev = bps.back();
	int i = bps.size() - 2; 
	while (i >= 0) {
		if (bps[i] + fopts.clusterTrimedge >= prev)
			trimInfo[i] = trimInfo[i + 1];
		else
			trimInfo[i] = i;
		prev = bps[i];
		i--;
	}
}

void modifyClusterBoundaries (vector<cluster> &clusts, uint32_t original, uint32_t after, bool axis) 
{
	if (after == original) return;
	uint32_t c;
	for (c = 0; c < clusts.size(); ++c) {
		assert(after > original);
		if (axis == 0) { // y-axis
			if (clusts[c].yStart == original) {
				clusts[c].yStart = after;
				clusts[c].xStart += (after - original);
			}
			if (clusts[c].yEnd == original) {
				clusts[c].yEnd = after;
				clusts[c].xEnd += (after - original);
			}			
		}
		else { // x-axis
			if (clusts[c].xStart == original) {
				clusts[c].xStart = after;
				clusts[c].yStart += (after - original);
			}
			if (clusts[c].xEnd == original) {
				clusts[c].xEnd = after;
				clusts[c].yEnd += (after - original);
			}
		}
	}
}

void trimClusters_nomodifybreakpoints (sample &a, vector<uint32_t> &bps, vector<uint32_t> &trimInfo, vector<cluster> &dense_clusts, 
					vector<cluster> &sparse_clusts, vector<bool> &remove, bool axis, bool unifysample = 0) 
{
	uint32_t i;
	for (i = 0; i < trimInfo.size(); ++i) {
		if (trimInfo[i] == i) continue;
		remove[i] = 1;
		modifyClusterBoundaries(dense_clusts, bps[i], bps[trimInfo[i]], axis);
		modifyClusterBoundaries(sparse_clusts, bps[i], bps[trimInfo[i]], axis);
		if (unifysample) {
			modifyClusterBoundaries(a.dense_clusts, bps[i], bps[trimInfo[i]], axis);
			modifyClusterBoundaries(a.sparse_clusts, bps[i], bps[trimInfo[i]], axis);			
		}
	}
}

void trimClusters (sample &a, vector<uint32_t> &bps, vector<uint32_t> &trimInfo, vector<cluster> &dense_clusts, 
					vector<cluster> &sparse_clusts, bool axis, bool unifysample = 0)
{
	vector<bool> remove(bps.size(), 0);
	trimClusters_nomodifybreakpoints (a, bps, trimInfo, dense_clusts, sparse_clusts, remove, axis, unifysample);
	REsize(bps, remove);
}

void collectBreakpoints (vector<cluster> &dense_clusts, vector<cluster> &sparse_clusts, set<uint32_t> & bpset, vector<uint32_t> &bps, bool axis)
{
	uint32_t j;
	for (j = 0; j < dense_clusts.size(); ++j) {
		if (axis == 0) {
			bpset.insert(dense_clusts[j].yStart);
			bpset.insert(dense_clusts[j].yEnd);			
		}
		else {
			bpset.insert(dense_clusts[j].xStart);
			bpset.insert(dense_clusts[j].xEnd);				
		}
	}
	for (j = 0; j < sparse_clusts.size(); ++j) {
		if (axis == 0) {
			bpset.insert(sparse_clusts[j].yStart);
			bpset.insert(sparse_clusts[j].yEnd);				
		}
		else {
			bpset.insert(sparse_clusts[j].xStart);
			bpset.insert(sparse_clusts[j].xEnd);	
		}
	}
  	for (auto it = bpset.begin(); it != bpset.end(); ++it)
    	bps.push_back(*it);
}

void firstTrim (sample &a, set<uint32_t> &bpset, vector<cluster> &dense_clusts, vector<cluster> &sparse_clusts, vector<uint32_t> &bps, 
				const fragopt_t &fopts, bool axis, bool unifysample = 0)
{
    // collect breakpoints from the axis (y-axis: 0; x-axis: 1)
    collectBreakpoints(dense_clusts, sparse_clusts, bpset, bps, axis);
    bpset.clear();

	// trim the breakpoints on the axis
	vector<uint32_t> trimInfo(bps.size());
	iota(trimInfo.begin(), trimInfo.end(), 0);
	clusterBreakpoints(bps, trimInfo, fopts);
	trimClusters(a, bps, trimInfo, dense_clusts, sparse_clusts, axis, unifysample);

	// if (fopts.debug and unifysample) {
	// 	ofstream fclust("trimlines_unify.bed", ios_base::app);
	//   	for (auto& it : bps)
	//     	fclust << it << "\t" << *(a.readname) << endl;	
	// 	fclust.close();			
	// }
}

// void secondTrim (sample &a, vector<cluster> &dense_clusts, vector<cluster> &sparse_clusts, vector<uint32_t> &bps_first, vector<uint32_t> &bps_second,
// 				 const fragopt_t &fopts, bool axis)
// {
// 	set<uint32_t> breakpoints_f;
// 	set<uint32_t> breakpoints_s;
// 	uint32_t i, diff;
//   	for (i = 0; i != bps_first.size(); ++i)
//     	breakpoints_f.insert(bps_first[i]);

// 	auto its = breakpoints_f.begin();
// 	auto ite = breakpoints_f.begin();
// 	for (i = 0; i < dense_clusts.size(); ++i) {
// 		if (axis == 0) {
// 			its = breakpoints_f.lower_bound(dense_clusts[i].yStart); 
// 			ite = breakpoints_f.lower_bound(dense_clusts[i].yEnd);
// 			assert(*its == dense_clusts[i].yStart and *ite == dense_clusts[i].yEnd);			
// 		}
// 		else {
// 			its = breakpoints_f.lower_bound(dense_clusts[i].xStart); 
// 			ite = breakpoints_f.lower_bound(dense_clusts[i].xEnd);
// 			assert(*its == dense_clusts[i].xStart and *ite == dense_clusts[i].xEnd);				
// 		}

// 		// project the first collected breakpoints to the other axis
// 		for (auto it = its; it != ite; ++it)
// 			diff = (*it > *its) ? (*it - *its) : 0;
// 			if (axis == 0) 
// 				breakpoints_s.insert(diff + dense_clusts[i].xStart);
// 			else
// 				breakpoints_s.insert(diff + dense_clusts[i].yStart);
// 	}

// 	for (i = 0; i < sparse_clusts.size(); ++i) {
// 		if (axis == 0) {
// 			its = breakpoints_f.lower_bound(sparse_clusts[i].yStart); 
// 			ite = breakpoints_f.lower_bound(sparse_clusts[i].yEnd);
// 			assert(*its == sparse_clusts[i].yStart and *ite == sparse_clusts[i].yEnd);			
// 		}
// 		else {
// 			its = breakpoints_f.lower_bound(sparse_clusts[i].xStart); 
// 			ite = breakpoints_f.lower_bound(sparse_clusts[i].xEnd);
// 			assert(*its == sparse_clusts[i].xStart and *ite == sparse_clusts[i].xEnd);				
// 		}

// 		// project the first collected breakpoints to the other axis		
// 		for (auto it = its; it != ite; ++it)
// 			diff = (*it > *its) ? (*it - *its) : 0;
// 			if (axis == 0)
// 				breakpoints_s.insert(diff + sparse_clusts[i].xStart);	
// 			else
// 				breakpoints_s.insert(diff + sparse_clusts[i].yStart);	
// 	}

//   	for (i = 0; i != bps_first.size(); ++i)
//     	breakpoints_s.insert(bps_first[i]);

//   	for (auto it = breakpoints_s.begin(); it != breakpoints_s.end(); ++it)
//     	bps_second.push_back(*it);

//     breakpoints_s.clear();
//     breakpoints_f.clear();

// 	// trim the breakpoints on the other axis (1 ^ axis)
// 	vector<uint32_t> trimInfo(bps_second.size());
// 	iota(trimInfo.begin(), trimInfo.end(), 0);
// 	clusterBreakpoints(bps_second, trimInfo, fopts);
// 	trimClusters(a, bps_second, trimInfo, dense_clusts, sparse_clusts, (1 ^ axis));
// 	// cerr << " get all the breakpoints on y, x-axis!" << endl;
// }

void project_helper (set<uint32_t>::iterator its, set<uint32_t>::iterator ite, cluster &clust, set<uint32_t> &bps, bool axis)
{
	// project the first collected breakpoints to the other axis
	uint32_t diff;	
			
	if (clust.strand == 0) {
		for (auto it = its; it != ite; ++it) {
			diff = (*it > *its) ? (*it - *its) : 0;
			// if (diff == 0)
			// 	continue;
			if (axis == 0) 
				bps.insert(diff + clust.xStart);
			else
				bps.insert(diff + clust.yStart);
		}		
	}
	else {
		for (auto it = ite; it != its; --it) {
			diff = (*ite > *it) ? (*ite - *it) : 0;
			// if (diff == 0)
			// 	continue;			
			if (axis == 0)
				bps.insert(diff + clust.xStart);
			else				
				bps.insert(diff + clust.yStart);
		}
	}
}


// project sample a to sample b
// axis == 0: project y-axis to x-axis
// axis == 1: project x-axis to y-axis
void project (sample &a, sample &b, vector<cluster> &dense_clusts, vector<cluster> &sparse_clusts, const fragopt_t &fopts, bool axis, bool dense)
{
	uint32_t i;
	set<uint32_t> bps;
	set<uint32_t> proj_bps;
	// insert original breakpoints of a to bps
  	for (auto&it : a.breakpoints) 
    	bps.insert(it);

	auto its = bps.begin();
	auto ite = bps.begin();

	// insert original clusters boundaries to bps
	if (dense) {
		for (auto&it : dense_clusts) {
			if (axis == 0) {
				bps.insert(it.yStart);
				bps.insert(it.yEnd);
			}
			else{
				bps.insert(it.xStart);
				bps.insert(it.xEnd);
			}
		}	

		for (i = 0; i < dense_clusts.size(); ++i) {
			if (axis == 0) {
				its = bps.lower_bound(dense_clusts[i].yStart); 
				ite = bps.lower_bound(dense_clusts[i].yEnd);
				assert(*its == dense_clusts[i].yStart and *ite == dense_clusts[i].yEnd);			
			}
			else {
				its = bps.lower_bound(dense_clusts[i].xStart); 
				ite = bps.lower_bound(dense_clusts[i].xEnd);
				assert(*its == dense_clusts[i].xStart and *ite == dense_clusts[i].xEnd);				
			}

			// project the first collected breakpoints to the other axis
			project_helper(its, ite, dense_clusts[i], proj_bps, axis);
		}		
	}
	else {
		for (auto&it : sparse_clusts) {
			if (axis == 0) {
				bps.insert(it.yStart);
				bps.insert(it.yEnd);
			}
			else{
				bps.insert(it.xStart);
				bps.insert(it.xEnd);
			}
		}	

		for (i = 0; i < sparse_clusts.size(); ++i) {
			if (axis == 0) {
				its = bps.lower_bound(sparse_clusts[i].yStart); 
				ite = bps.lower_bound(sparse_clusts[i].yEnd);
				assert(*its == sparse_clusts[i].yStart and *ite == sparse_clusts[i].yEnd);			
			}
			else {
				its = bps.lower_bound(sparse_clusts[i].xStart); 
				ite = bps.lower_bound(sparse_clusts[i].xEnd);
				assert(*its == sparse_clusts[i].xStart and *ite == sparse_clusts[i].xEnd);				
			}

			// project the first collected breakpoints to the other axis	
			project_helper(its, ite, sparse_clusts[i], proj_bps, axis);
		}	
	}

	// set<uint32_t> b_sort;
	for (auto&it : b.breakpoints)
		proj_bps.insert(it);

	b.breakpoints.clear();
  	for (auto it = proj_bps.begin(); it != proj_bps.end(); ++it)
    	b.breakpoints.push_back(*it);	
}


// project sample a to sample b
// axis == 0: project y-axis to x-axis
// axis == 1: project x-axis to y-axis
void selfproject (sample &a, vector<cluster> &dense_clusts, vector<cluster> &sparse_clusts, const fragopt_t &fopts, bool dense)
{
	uint32_t i;
	set<uint32_t> bps;
	set<uint32_t> temp_bps;
	// insert original breakpoints of a to bps
  	for (auto&it : a.breakpoints) 
    	bps.insert(it);

	auto its = bps.begin();
	auto ite = bps.begin();

	// insert original clusters boundaries to bps
	if (dense) {
		for (auto&it : dense_clusts) {
			bps.insert(it.yStart);
			bps.insert(it.yEnd);
		}

		for (i = 0; i < dense_clusts.size(); ++i) {
			its = bps.lower_bound(dense_clusts[i].yStart); 
			ite = bps.lower_bound(dense_clusts[i].yEnd);
			assert(*its == dense_clusts[i].yStart and *ite == dense_clusts[i].yEnd);			

			// project the first collected breakpoints to the other axis
			temp_bps.clear();
			project_helper(its, ite, dense_clusts[i], temp_bps, 0);

			for (auto&tmp : temp_bps) 
				bps.insert(tmp);
		}		
	}
	else {
		for (auto&it : sparse_clusts) {
			bps.insert(it.yStart);
			bps.insert(it.yEnd);
		}	

		for (i = 0; i < sparse_clusts.size(); ++i) {
			its = bps.lower_bound(sparse_clusts[i].yStart); 
			ite = bps.lower_bound(sparse_clusts[i].yEnd);
			assert(*its == sparse_clusts[i].yStart and *ite == sparse_clusts[i].yEnd);			

			// project the first collected breakpoints to the other axis	
			temp_bps.clear();
			project_helper(its, ite, sparse_clusts[i], temp_bps, 0);
			for (auto&tmp : temp_bps) 
				bps.insert(tmp);
		}	
	}

	a.breakpoints.clear();
  	for (auto it = bps.begin(); it != bps.end(); ++it)
    	a.breakpoints.push_back(*it);	
}

void secondTrim (sample &a, vector<cluster> &dense_clusts, vector<cluster> &sparse_clusts, vector<uint32_t> &bps_first, vector<uint32_t> &bps_second,
				 const fragopt_t &fopts, bool axis, bool unifysample)
{
	set<uint32_t> breakpoints_f;
	set<uint32_t> breakpoints_s;
	uint32_t i;
  	for (i = 0; i != bps_first.size(); ++i) {
    	breakpoints_f.insert(bps_first[i]);
  	}

	auto its = breakpoints_f.begin();
	auto ite = breakpoints_f.begin();
	for (i = 0; i < dense_clusts.size(); ++i) {
		if (axis == 0) {
			its = breakpoints_f.lower_bound(dense_clusts[i].yStart); 
			ite = breakpoints_f.lower_bound(dense_clusts[i].yEnd);
			assert(*its == dense_clusts[i].yStart and *ite == dense_clusts[i].yEnd);			
		}
		else {
			its = breakpoints_f.lower_bound(dense_clusts[i].xStart); 
			ite = breakpoints_f.lower_bound(dense_clusts[i].xEnd);
			assert(*its == dense_clusts[i].xStart and *ite == dense_clusts[i].xEnd);				
		}

		// project the first collected breakpoints to the other axis
		project_helper(its, ite, dense_clusts[i], breakpoints_s, axis);
	}

	for (i = 0; i < sparse_clusts.size(); ++i) {
		if (axis == 0) {
			its = breakpoints_f.lower_bound(sparse_clusts[i].yStart); 
			ite = breakpoints_f.lower_bound(sparse_clusts[i].yEnd);
			assert(*its == sparse_clusts[i].yStart and *ite == sparse_clusts[i].yEnd);			
		}
		else {
			its = breakpoints_f.lower_bound(sparse_clusts[i].xStart); 
			ite = breakpoints_f.lower_bound(sparse_clusts[i].xEnd);
			assert(*its == sparse_clusts[i].xStart and *ite == sparse_clusts[i].xEnd);				
		}

		// project the first collected breakpoints to the other axis	
		project_helper(its, ite, sparse_clusts[i], breakpoints_s, axis);
	}

	if (unifysample == 0) {
	  	for (i = 0; i != bps_first.size(); ++i)
	    	breakpoints_s.insert(bps_first[i]);		
	}

  	for (auto it = breakpoints_s.begin(); it != breakpoints_s.end(); ++it)
    	bps_second.push_back(*it);

    breakpoints_s.clear();

	if (unifysample == 1)
		return;

    bps_first.clear();

	// trim the breakpoints on the other axis (1 ^ axis)
	vector<uint32_t> trimInfo(bps_second.size());
	iota(trimInfo.begin(), trimInfo.end(), 0);
	clusterBreakpoints(bps_second, trimInfo, fopts);
	trimClusters(a, bps_second, trimInfo, dense_clusts, sparse_clusts, (1 ^ axis));
	// cerr << " get all the breakpoints on y, x-axis!" << endl;
	return;
}

// cluster breakpoints in edge around pivot
// No change to the actual boundaries of clusters
// return the updated breakpoints for sample a 
void updateBreakpointsBasedOnPivot(vector<uint32_t> &breakpoints, vector<uint32_t> &pivot, vector<uint32_t> &edge, const fragopt_t &fopts) 
{
	set<uint32_t> p;
	set<uint32_t> e;	
	uint32_t i;
  	for (i = 0; i != pivot.size(); ++i)
    	p.insert(pivot[i]);

  	for (i = 0; i != edge.size(); ++i) 
    	e.insert(edge[i]);

    edge.clear();  

	for (auto it = e.begin(); it != e.end(); ++it) 		
		edge.push_back(*it);

	vector<uint32_t> trimInfo(edge.size());
	iota(trimInfo.begin(), trimInfo.end(), 0);

	auto low = p.begin();
	auto high = p.end();
	uint32_t l, h, temp;
	uint32_t m = numeric_limits<uint32_t>::max();
	for (i = 0; i < edge.size(); ++i) {
		l = edge[i] > fopts.clusterTrimedge ? edge[i] - fopts.clusterTrimedge : 0;
		h = edge[i] + fopts.clusterTrimedge;
		low = p.lower_bound(l);
		high = p.lower_bound(h);
		for (auto it = low; it != high; ++it) {
			temp = *it > edge[i] ? *it - edge[i] : edge[i] - *it;
			m = min(m, temp);
			if (temp == m) 				
				trimInfo[i] = distance(it, p.begin());
		}
	}

	// update sample a's breakpoints
	for (i = 0; i < trimInfo.size(); ++i) {
		if (trimInfo[i] = i)
			p.insert(edge[i]);
	}		
	breakpoints.clear();
	for (auto it = p.begin(); it != p.end(); ++it) 
		breakpoints.push_back(*it);
}









