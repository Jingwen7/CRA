#ifndef BREAKPOINT_H_
#define BREAKPOINT_H_

#include <numeric>
#include <vector>
#include "cluster.h"

void trimBreakpoints(const vector<uint32_t> &bps, vector<uint32_t> &trimInfo);

void modifyClusterBoundaries (vector<cluster> &clusts, uint32_t original, uint32_t after);

void trimClusters(string * readname, vector<uint32_t> &bps, vector<uint32_t> &trimInfo, vector<cluster> &dense_clusts, vector<cluster> &sparse_clusts);

void trimOnY(string * readname, vector<cluster> &dense_clusts, vector<cluster> &sparse_clusts, vector<uint32_t> &bps, const fragopt_t &fopts);

void trimOnX(string * readname, vector<cluster> &dense_clusts, vector<cluster> &sparse_clusts, vector<uint32_t> &bps_Y, vector<uint32_t> &bps_X, const fragopt_t &fopts);

#endif