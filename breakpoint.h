#ifndef BREAKPOINT_H_
#define BREAKPOINT_H_

#include <numeric>
#include <vector>
#include "cluster.h"
#include "sample.h"

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

#endif