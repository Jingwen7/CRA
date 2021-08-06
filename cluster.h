#ifndef CLUSTER_H_
#define CLUSTER_H_

#include "rfpriv.h"
#include "option.h"
#include "hit.h"

class cluster
{
public:
	bool strand;
	uint32_t s, e;
	uint32_t xStart, xEnd, yStart, yEnd;
	uint32_t cidx;

	cluster (uint32_t S, uint32_t E, uint32_t XStart, uint32_t XEnd, uint32_t YStart, uint32_t YEnd, bool st) : 
	s(S), e(E), xStart(XStart), xEnd(XEnd), yStart(YStart), yEnd(YEnd), strand(st) {};

	cluster (uint32_t XStart, uint32_t XEnd, uint32_t YStart, uint32_t YEnd, bool st, uint32_t M) : 
	xStart(XStart), xEnd(XEnd), yStart(YStart), yEnd(YEnd), strand(st), cidx(M) {};
	~cluster () {};
};

bool clustDiagonalSort (const cluster &a, const cluster &b); 

bool clustAntiDiagonalSort (const cluster &a, const cluster &b);

class IntervalSet 
{
public:
	double slope;
	double intercept;
	bool strand;
	vector<pair<uint32_t, bool>> Set;

	IntervalSet (cluster & clust) {
		slope = (double)((int64_t)clust.yEnd - (int64_t)clust.yStart)/((int64_t)clust.xEnd - (int64_t)clust.xStart);
		if (clust.strand == 0) {
			intercept = ((double)((int64_t)clust.xEnd * clust.yStart - (int64_t)clust.xStart * clust.yEnd))/((int64_t) clust.xEnd - (int64_t) clust.xStart);
		}
		else {
			slope = -1 * slope;
			intercept = (double)((int64_t)clust.xStart * clust.yStart - (int64_t)clust.xEnd * clust.yEnd)/((int64_t)clust.xStart - (int64_t)clust.xEnd);
		}
		strand = clust.strand;
	};

	~IntervalSet() {};

	int operator() (const pair<uint32_t, bool> &a, const pair<uint32_t, bool> &b) {
		if (a.second == b.second and a.second == 0) {
			return a.first < b.first;
		}
		else if (a.second == b.second and a.second == 1) {
			if (strand == 0) return a.first < b.first;
			else return a.first > b.first;
		}
		else if (a.second == 0 and b.second == 1) {
			if (strand == 0) return a.first*slope + intercept < (double) b.first;
			else return a.first*slope + intercept > (double) b.first;
		}
		else {
			if (strand == 0) return (double) a.first < b.first*slope + intercept;
			else return (double) a.first > b.first*slope + intercept;
		}		
	}

	void Sort () {
		sort(Set.begin(), Set.end(), *this);
	}
};

typedef vector<cluster> clusters;

void storeDiag(vector<hit> &hits, clusters &clust, idxopt_t &iopt, fragopt_t &fopt, bool st = 0);

void storeDiagCluster(vector<hit> &hits, clusters &clust, bool st, idxopt_t &iopt, fragopt_t &fopt);

void mergeDiagCluster(clusters &clusts, int64_t mergeDiag, bool st);

void trim(clusters &clusts, const set<pointInfo> &Spoints, const set<pointInfo> &Epoints, int clusterTrimedge);

void trim_ovpClusters(clusters &clusts, int clusterTrimedge);

void splitClusters(clusters & clusts, clusters & splitclusts);

void fragLabel(const clusters &clusts, fragopt_t &fopts);




#endif