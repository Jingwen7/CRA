#ifndef HIT_H
#define HIT_H

#include "rfpriv.h"
#include "cluster.h"
#include "index.h"
#include "sample.h"
#include <vector>

class hit {
public:
	uint256_t mi;
	uint32_t x; 
	uint32_t y;

	hit() {};
	hit(uint32_t &Mi, uint32_t X, uint32_t Y) : mi(Mi), x(X), y(Y) {};
	~hit() {};	
};

void rf_hit(const sample &sample_a, const sample &sample_b, const uint32_t a, const uint32_t b, const idx_t &mi, vector<hit> &fhits, 
				vector<hit> &rhits, const fragopt_t &fopts, const idxopt_t &iopt, bool self);

void storeDiagCluster (uint32_t s, uint32_t e, vector<hit> &hits, clusters &clust, bool st, const idxopt_t &iopt, const fragopt_t &fopts);

void cleanDiag(vector<hit> &hits, clusters &clust, const fragopt_t &fopts, const idxopt_t &iopt, bool st);

inline int64_t DiagonalDifference(hit &a, hit &b, bool strand = 0) {
	int64_t aDiag, bDiag;
	if (strand == 0) { // Matches
		aDiag = (int64_t)a.x - (int64_t)a.y;
		bDiag = (int64_t)b.x - (int64_t)b.y;
		return aDiag - bDiag;		
	}
	else { // revMathches
		aDiag = (int64_t) a.x + (int64_t) a.y; 
		bDiag= (int64_t) b.x + (int64_t) b.y;
		return aDiag - bDiag;				
	}
}


inline int64_t maxGapDifference(hit &a, hit &b) {
	int64_t aDiff = abs((int64_t)b.x - (int64_t)a.x);
	int64_t bDiff = abs((int64_t)b.y - (int64_t)a.y);
	return max(aDiff, bDiff);
}

inline int64_t minGapDifference (hit &a, hit &b) {
	int64_t aDiff = abs((int64_t)b.x - (int64_t)a.x);
	int64_t bDiff = abs((int64_t)b.y - (int64_t)a.y);
	return min(aDiff, bDiff);
}

bool hitDiagonalSort (const hit &a, const hit &b);

bool hitAntiDiagonalSort (const hit &a, const hit &b);	

#endif