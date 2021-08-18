#ifndef BREAKPOINT_H_
#define BREAKPOINT_H_

#include <numeric>
#include <vector>
#include "cluster.h"
#include "sample.h"

void checkBreakpoints_Clusters (vector<uint32_t> &bps, vector<cluster> &clusts, bool axis);

void checkBreakpoints_Clusters (set<uint32_t> &s, vector<cluster> &clusts, bool axis);

void insertPoint (vector<cluster> &clusts, set<uint32_t> &s, bool axis);

void insertPoint (vector<cluster> &clusts, vector<uint32_t> &s, bool axis);

template<typename T>
void REsize(vector<T> &original, const vector<bool> &remove);


void clusterBreakpoints (const vector<uint32_t> &bps, vector<uint32_t> &trimInfo, const fragopt_t &fopts);

void modifyClusterBoundaries (vector<cluster> &clusts, uint32_t original, uint32_t after, bool axis);

void collectBreakpoints (vector<cluster> &dense_clusts, vector<cluster> &sparse_clusts, vector<uint32_t> &bps, bool axis);

void trimClusters_nomodifybreakpoints (sample &a, vector<uint32_t> &bps, vector<uint32_t> &trimInfo, vector<cluster> &dense_clusts, 
					vector<cluster> &sparse_clusts, vector<bool> &remove, bool axis, bool unifysample);

void trimClusters (sample &a, vector<uint32_t> &bps, vector<uint32_t> &trimInfo, vector<cluster> &dense_clusts, 
					vector<cluster> &sparse_clusts, bool axis, bool unifysample);
/*
This function collect breakpoints (on one axis) from the clusters and clustering the breakpoints;
meanwhile make necessary trimming to the clusters
*/
void firstTrim (sample &a, set<uint32_t> &bpset, vector<cluster> &dense_clusts, vector<cluster> &sparse_clusts, vector<uint32_t> &bps, 
				const fragopt_t &fopts, bool axis, bool unifysample);
/*
This function projected the previously collected breakpoints (from one axis) "bps_first" to the other axis "bps_second";
clustering the breakpoints on the other axis "bps_second";
meanwhile make necessary trimming to the clusters;
*/
void secondTrim (sample &a, vector<cluster> &dense_clusts, vector<cluster> &sparse_clusts, vector<uint32_t> &bps_first, 
				vector<uint32_t> &bps_second, const fragopt_t &fopts, bool axis, bool unifysample);

void updateBreakpointsBasedOnPivot(sample &a, vector<uint32_t> &pivot, vector<uint32_t> &edge, const fragopt_t &fopts);

void project (sample &a, sample &b, vector<cluster> &dense_clusts, vector<cluster> &sparse_clusts, const fragopt_t &fopts, bool axis);

#endif