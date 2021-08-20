#include "rfpriv.h"
#include "hit.h"
#include "index.h"
#include "cluster.h"
#include <vector>
#include <map>

static uint256_t rmask = (uint256_t) 1; // 00000...1
static uint256_t fmask = (((uint256_t)1 << 255) - 1) << 1; // 111...1110

bool hitDiagonalSort (const hit &a, const hit &b)
{
	int64_t aDiag = (int64_t) a.x - (int64_t) a.y;
	int64_t bDiag = (int64_t) b.x - (int64_t) b.y;
	if (aDiag != bDiag)
		return aDiag < bDiag;
	else
		return a.x < b.x;
}	

bool hitAntiDiagonalSort (const hit &a, const hit &b)
{
	int64_t aDiag = (int64_t) a.x + (int64_t) a.y;
	int64_t bDiag = (int64_t) b.x + (int64_t) b.y;
	if (aDiag != bDiag)
		return aDiag < bDiag;
	else
		return a.x < b.x;
}	

bool hitCartesianSort (const hit &a, const hit &b)
{
	if (a.x != b.x)
		return a.x < b.x;
	else
		return a.y < b.y;
}	

bool hitAntiCartesianSort (const hit &a, const hit &b)
{
	if (a.x != b.x)
		return a.x < b.x;
	else
		return a.y > b.y;
}

inline void find_match (vector<hit> &hits, const uint256_t &mi, const vector<uint32_t> &pos, const kmerhash &khash, bool self)
{
	uint32_t i, j;
	hit temp;
	if (self) { // only half of the self-self alignment
		for (i = 0; i < pos.size(); ++i) {
			for (j = i + 1; j < pos.size(); ++j) {
				temp.mi = (mi & fmask); 
				temp.x = pos[i]; 
				temp.y = pos[j];
				hits.push_back(temp);
			}
		}
	}
	else {
		auto ot = khash.find(mi);
		if (ot != khash.end()) {
			for (i = 0; i < ot->second.size(); ++i) {
				for (j = 0; j < pos.size(); ++j) {
					temp.mi = (mi & fmask); 
					temp.x = pos[j]; 
					temp.y = (ot->second)[i];
					hits.push_back(temp);	
				}
			}
		}		
	}
}

inline void find_match_boundary_checking (const sample &sample_i, const sample &sample_j, vector<hit> &hits, const uint256_t &mi, const vector<uint32_t> &pos, const kmerhash &khash, bool self)
{
	uint32_t i, j;
	hit temp;
	if (self) { // only half of the self-self alignment
		for (i = 0; i < pos.size(); ++i) {
			for (j = i + 1; j < pos.size(); ++j) {
				// if (pos[i] >= sample_i.breakpoints[0] and pos[i] <= sample_i.breakpoints.back()
					// and pos[j] >= sample_j.breakpoints[0] and pos[j] <= sample_j.breakpoints.back()) {
					temp.mi = (mi & fmask); 
					temp.x = pos[i]; 
					temp.y = pos[j];
					hits.push_back(temp);
				// }
			}
		}
	}
	else {
		auto ot = khash.find(mi);
		if (ot != khash.end()) {
			for (j = 0; j < ot->second.size(); ++j) {
				for (i = 0; i < pos.size(); ++i) {
					if (pos[i] >= sample_i.breakpoints[0] and pos[i] <= sample_i.breakpoints.back()
					and (ot->second)[j] >= sample_j.breakpoints[0] and (ot->second)[j] <= sample_j.breakpoints.back()) {
						temp.mi = (mi & fmask); 
						temp.x = pos[i]; 
						temp.y = (ot->second)[j];
						hits.push_back(temp);	
					}
				}
			}
		}		
	}
}

void checkForwardmatch (const Genome *genome, uint32_t a, uint32_t b, uint32_t aPos, uint32_t bPos, const idxopt_t &iopt)
{
	if (strncmp(&(genome->seqs[a])[aPos], &(genome->seqs[b])[bPos], iopt.k) != 0) {
		cerr << "match is wrong" << endl;
	}
}

void rf_hit(const vector<sample> & samples, const uint32_t a, const uint32_t b, const idx_t &mi, vector<hit> &fhits, vector<hit> &rhits, const fragopt_t &fopts, const idxopt_t &iopt, bool self = 1) 
{
	// cerr << "mi.seqh[a].size(): " <<  mi.seqh[a]->size() << endl;
	uint256_t f, r; bool z;
	for (auto it = mi.seqh[a]->begin(); it != mi.seqh[a]->end(); ++it) {
		// minimizer: it->first;
		// vector<pos>: it->second; 
		if (it->second.size() < fopts.freq) {
			z = ((rmask & it->first) == 1) ? 1 : 0; // direction
			f = (it->first & fmask);
			r = (it->first | rmask);
			if (self) { // self-self only has forward matches
				if ((z ^ 0) == 0) // z = 0: forward 
					find_match(fhits, f, it->second, *mi.seqh[b], self);
				else if ((z ^ 1) == 0) // z = 1: forward
					find_match(fhits, r, it->second, *mi.seqh[b], self);
			}
			else { // sample i to sample j matches
				if ((z ^ 0) == 0) // z = 0: forward 
					find_match_boundary_checking(samples[a], samples[b], fhits, f, it->second, *mi.seqh[b], self);
				if ((z ^ 1) == 0) // z = 1: forward
					find_match_boundary_checking(samples[a], samples[b], fhits, r, it->second, *mi.seqh[b], self);
				if ((z ^ 1) == 1) // z = 0: reverse matches
					find_match_boundary_checking(samples[a], samples[b], rhits, r, it->second, *mi.seqh[b], self);
				if ((z ^ 0) == 1) // z = 1: reverse matches
					find_match_boundary_checking(samples[a], samples[b], rhits, f, it->second, *mi.seqh[b], self);
			}
			
		}
	}
	// if (fopts.debug and !self) {
	// 	ofstream fclust("hit.bed", ios::app);
	// 	cerr << "forward hits: " << fhits.size() << endl;
	// 	for (uint32_t m = 0; m < fhits.size(); ++m) {
	// 		// checkForwardmatch(mi.genome, a, b, fhits[m].x, fhits[m].y, iopt);
	// 		fclust << fhits[m].x << "\t" << fhits[m].y << "\t" << fhits[m].x + iopt.k << "\t" << fhits[m].y + iopt.k << "\t" << iopt.k << "\t" << "0" << endl;
	// 	}
	// 	cerr << "reverse hits: " << rhits.size() << endl;
	// 	for (uint32_t m = 0; m < rhits.size(); ++m) {
	// 		fclust << rhits[m].x << "\t" << rhits[m].y << "\t" << rhits[m].x + iopt.k << "\t" << rhits[m].y + iopt.k << "\t" << iopt.k << "\t" << "1" << endl;
	// 	}
	// 	fclust.close();			
	// }
}

bool CloseToPreviousCluster(const cluster &prev, uint32_t xStart, uint32_t xEnd, uint32_t yStart, uint32_t yEnd, const fragopt_t &fopts) {
	int64_t xDiff = abs((int64_t)xStart - (int64_t)prev.xEnd);
	int64_t yDiff;
	if (prev.strand == 0) 
		yDiff = abs((int64_t) yStart - (int64_t) prev.yEnd);
	else 
		yDiff = abs((int64_t) prev.yStart - (int64_t) yEnd);	// (TODO): make sure this is working for reversed clusters!

	int64_t curDiag = 0, prevDiag;
	if (prev.strand == 0) {
		curDiag = (int64_t) yStart - (int64_t) xStart;
		prevDiag = (int64_t) prev.yEnd - (int64_t) prev.xEnd;
	}
	else {
		curDiag = (int64_t) xStart + (int64_t) yEnd;
		prevDiag = (int64_t) prev.xStart + (int64_t) prev.yEnd;
	}

	if (max(xDiff, yDiff) <= 100 and abs(curDiag - prevDiag) <= 50) 
		return true;
	return false;
}

void UpdateCluster(cluster &prev, uint32_t xS, uint32_t xE, uint32_t yS, uint32_t yE) {
	prev.xStart = min(prev.xStart, xS);
	prev.xEnd = max(prev.xEnd, xE);
	prev.yStart = min(prev.yStart, yS);
	prev.yEnd = max(prev.yEnd, yE);
}

void MergeTwoClusters(cluster &prev, uint32_t xS, uint32_t xE, uint32_t yS, uint32_t yE, uint32_t end) {
	UpdateCluster(prev, xS, xE, yS, yE);
	prev.e = end;	
}

void storeDiagCluster (uint32_t s, uint32_t e, vector<hit> &hits, clusters &clust, bool st, const idxopt_t &iopt, const fragopt_t &fopts) 
{
	sort(hits.begin() + s, hits.begin() + e, hitCartesianSort);
	uint32_t initial = clust.size();
	uint32_t cs = s, ce = s;
	int64_t diagdiff, gapdiff;
	while (cs < e) {
		ce = cs + 1;
		uint32_t  xStart = hits[cs].x, 
				  xEnd = hits[cs].x + iopt.k, 
			      yStart = hits[cs].y, 
			      yEnd = hits[cs].y + iopt.k;

		while (ce < e and abs(DiagonalDifference(hits[ce], hits[ce - 1], st)) < fopts.clusterMaxDiag and minGapDifference(hits[ce], hits[ce - 1]) <= fopts.clusterMaxDist ) {
			xStart = min(xStart, hits[ce].x);
			xEnd   = max(xEnd, hits[ce].x + iopt.k);
			yStart = min(yStart, hits[ce].y);
			yEnd   = max(yEnd, hits[ce].y + iopt.k);
			ce++;
		}	

		if (ce - cs >= fopts.DiagMinCluster and xEnd - xStart >= fopts.clusterMinLength and yEnd - yStart >= fopts.clusterMinLength and (float) (ce - cs) / (xEnd - xStart) >= 0.05) {
			if (clust.size() > 0 and CloseToPreviousCluster(clust.back(), xStart, xEnd, yStart, yEnd, fopts)) {
				MergeTwoClusters(clust.back(), xStart, xEnd, yStart, yEnd, ce);
			}
			else {
				clust.push_back(cluster(cs, ce, xStart, xEnd, yStart, yEnd, st));
			}
		}
		cs = ce;
	}	

	// if (fopts.debug) {
	// 	ofstream rclust("diagcluster.bed", ios_base::app);
	// 	uint32_t m;
	// 	for (m = initial; m < clust.size(); m++) {
	// 		for (int c = clust[m].s; c < clust[m].e; ++c) {
	// 			rclust << hits[c].x << "\t" << hits[c].y << "\t" << hits[c].x + iopt.k << "\t"
	// 				   << hits[c].y + iopt.k << "\t" << iopt.k << "\t"  << clust[m].strand  << "\t" << m << endl;				
	// 		}
	// 	}
	// 	rclust.close();		
	// }
}

void cleanDiag(vector<hit> &hits, clusters &clust, const fragopt_t &fopts, const idxopt_t &iopt, bool st)
{
	//sort hits
	if (!st) {
		sort(hits.begin(), hits.end(), hitDiagonalSort);
	}
	else {
		sort(hits.begin(), hits.end(), hitAntiDiagonalSort);
	}
	uint32_t n = hits.size();
	vector<bool> onDiag(n, 0);
	vector<int> counts(n, -1);
	vector<int> rev_counts(n, -1);
	int counter = 0;
	if (n <= 1) return;
	//
	// starting from forward order
	//
	vector<bool> foward_onDiag(n, false);
	vector<bool> reverse_onDiag(n, false);

	uint32_t i, j;
	for (i = 1; i < n; ++i) {
		if (abs(DiagonalDifference(hits[i], hits[i - 1], st)) < fopts.CleanMaxDiag)	
			foward_onDiag[i - 1] = true;
	}

	bool prevOnDiag = false; int diagStart;
	for (i = 0; i < n; ++i) {
		if (prevOnDiag == false and foward_onDiag[i] == true) {
			diagStart = i;
			prevOnDiag = foward_onDiag[i];
		}
		else if (prevOnDiag == true and foward_onDiag[i] == false) {
			prevOnDiag = false;
			if (i - diagStart + 1 < fopts.DiagMinCluster)// [diagStart, i]
				for (j = diagStart; j <= i; ++j) { foward_onDiag[j] = false; }
			else {
				foward_onDiag[i] = true; 
				for (j = diagStart; j <= i; ++j) { counts[j] = counter; }
				counter++;
			}
		}
		else
			prevOnDiag = foward_onDiag[i];
	}	

	//
	// starting from reverse order
	//
	for (i = n - 2; i >= 0; --i) {
		if (abs(DiagonalDifference(hits[i], hits[i + 1], st)) < fopts.CleanMaxDiag)
			reverse_onDiag[i + 1] = true;
		if (i == 0) break;
	}

	prevOnDiag = false;
	for (i = n - 1; i >= 0; --i) {
		if (prevOnDiag == false and reverse_onDiag[i] == true) {
			diagStart = i;
			prevOnDiag = reverse_onDiag[i];
		}
		else if (prevOnDiag == true and reverse_onDiag[i] == false) {
			prevOnDiag = false;
			if (diagStart - i + 1 < fopts.DiagMinCluster) {// [diagStart, i]
				for (j = i; j <= diagStart; ++j) { 
					reverse_onDiag[j] = false; 
					counts[j] = -1;
				}
			}
			else {
				reverse_onDiag[i] = true;
				counter++;
			}
		}
		else
			prevOnDiag = reverse_onDiag[i];
		if (i == 0) break;
	}

	for (i = 0; i < n; ++i) {
		if (foward_onDiag[i] == true and reverse_onDiag[i] == true) 
			onDiag[i] = true; 
		else 
			onDiag[i] = false;
	}	

	uint32_t c = 0;
	for (i = 0; i < n; ++i) {
		if (onDiag[i]) {
			hits[c] = hits[i];
			counts[c] = counts[i];
			c++;
		}
	}
	cerr << "before cleaning: " << hits.size() << endl;
	cerr << "after cleaning: " << c << endl;
	hits.resize(c);
	counts.resize(c);

	// if (fopts.debug) {
	// 	ofstream fclust("cleanmatches.bed", ios_base::app);
	// 	for (uint32_t m = 0; m < hits.size(); ++m) {
	// 		// checkForwardmatch(mi.genome, a, b, fhits[m].x, fhits[m].y, iopt);
	// 		fclust << hits[m].x << "\t" << hits[m].y << "\t" << hits[m].x + iopt.k << "\t" << hits[m].y + iopt.k << "\t" << iopt.k << "\t" << st << "\t" << "999"<< endl;
	// 	}
	// 	fclust.close();			
	// }

	//
	// Store diagonal in clusters
	//
	uint32_t count_s = 0; c = 1;
	while (c <= counts.size()) {
		if (counts[c] == counts[c-1]) {c++; continue;}
		if (c == counts.size() and count_s == c) break;
		assert(count_s < c);
		storeDiagCluster (count_s, c, hits, clust, st, iopt, fopts);
		count_s = c;	
		c++;			
	}	
	if (count_s < c and c == counts.size() + 1) {
		storeDiagCluster (count_s, c - 1, hits, clust, st, iopt, fopts);		
	}
};
















