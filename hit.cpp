#include "rfpriv.h"
#include "hit.h"
#include "index.h"
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

inline void find_match (vector<hit> &hits, const uint256_t &mi, const vector<uint32_t> &pos, const kmerhash &khash, bool acrossStrand)
{
	uint32_t i, j;
	hit temp;
	if (!acrossStrand) { // only half of the self-self alignment
		for (i = 0; i < pos.size(); ++i) {
			for (j = i + 1; j < pos.size(); ++j) {
				// cerr << "i: " << i << " j: " << j << endl;
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

void checkForwardmatch (const Genome *genome, uint32_t a, uint32_t b, uint32_t aPos, uint32_t bPos, const idxopt_t &iopt)
{
	if (strncmp(&(genome->seqs[a])[aPos], &(genome->seqs[b])[bPos], iopt.k) != 0) {
		cerr << "match is wrong" << endl;
	}
}

void rf_hit(const uint32_t a, const uint32_t b, const idx_t &mi, vector<hit> &fhits, vector<hit> &rhits, const fragopt_t &fopts, const idxopt_t &iopt, bool self=1) 
{
	cerr << "mi.seqh[a].size(): " <<  mi.seqh[a]->size() << endl;
	// uint32_t count = 0;
	for (auto it = mi.seqh[a]->begin(); it != mi.seqh[a]->end(); ++it) {
		// minimizer: it->first;
		// vector<pos>: it->second; 
		// cerr << "minimizer: " << it->first << endl;
		// count++;
		uint256_t z, f, r;
		// uint256_t rmask, fmask, 
		// rmask = (uint256_t) 1; // 00000...1
		// fmask = (((uint256_t)1 << 255) - 1) << 1; // 111...1110
		z = ((rmask & it->first) == 1) ? 1 : 0; // direction
		f = (it->first & fmask);
		r = (it->first | rmask);
		// cerr << "it->second.size(): " << it->second.size() << endl;
		if (it->second.size() < fopts.freq) {
			if (z ^ 0 == 0) {// z = 0: forward 
				// cerr << "forward" << endl;
				find_match(fhits, f, it->second, *mi.seqh[b], 0);
			}
			else if (z ^ 1 == 0) {// z = 1: forward
				// cerr << "forward" << endl;			
				find_match(fhits, r, it->second, *mi.seqh[b], 0);
			}	
			else if (z ^ 1 == 1) {// z = 0: reverse
				// cerr << "reverse" << endl;
				find_match(rhits, r, it->second, *mi.seqh[b], 1);
			}
			else if (z ^ 0 == 1) {// reverse
				// cerr << "reverse" << endl;
				find_match(rhits, f, it->second, *mi.seqh[b], 1);
			}			
		}
	}
	if (fopts.debug) {
		ofstream fclust("matches.bed");
		cerr << "forward matches: " << fhits.size() << endl;
		for (uint32_t m = 0; m < fhits.size(); ++m) {
			// checkForwardmatch(mi.genome, a, b, fhits[m].x, fhits[m].y, iopt);
			fclust << fhits[m].x << "\t" << fhits[m].y << "\t" << fhits[m].x + iopt.k << "\t" << fhits[m].y + iopt.k << "\t" << "0" << "\t" << fhits[m].mi <<  endl;
		}
		cerr << "reverse matches: " << rhits.size() << endl;
		for (uint32_t m = 0; m < rhits.size(); ++m) {
			fclust << rhits[m].x << "\t" << rhits[m].y << "\t" << rhits[m].x + iopt.k << "\t" << rhits[m].y + iopt.k << "\t" << "1" << "\t" << rhits[m].mi << endl;
		}
		fclust.close();			
	}
}

void cleanDiag(vector<hit> &hits, const fragopt_t &fopts, bool st)
{
	uint32_t n = hits.size();
	vector<bool> onDiag(n, 0);
	if (n <= 1) return;
	//
	// starting from forward order
	//
	uint32_t nOnDiag2 = 0;
	vector<bool> foward_onDiag(n, false);
	vector<bool> reverse_onDiag(n, false);

	uint32_t i, j;
	for (i = 1; i < n; ++i) {
		if (abs(DiagonalDifference(hits[i], hits[i - 1], st)) < fopts.CleanMaxDiag) {	
			foward_onDiag[i - 1] = true;
			nOnDiag2++;
		}
	}

	bool prevOnDiag = false; int diagStart; int counter = 0;
	for (i = 0; i < n; ++i) {
		if (prevOnDiag == false and foward_onDiag[i] == true) 
			diagStart = i;
		if (prevOnDiag == true and foward_onDiag[i] == false) {
			if (i - diagStart + 1 < fopts.DiagMinCluster) { // [diagStart, i]
				for (j = diagStart; j <= i; ++j) { foward_onDiag[j] = false; }
			}
			else { foward_onDiag[i] = true; }
			counter++;
		}
		prevOnDiag = foward_onDiag[i];
	}	

	//
	// starting from reverse order
	//
	for (i = n - 2; i >= 0; --i) {
		if (abs(DiagonalDifference(hits[i], hits[i + 1], st)) < fopts.CleanMaxDiag) {	
			reverse_onDiag[i + 1] = true;
			nOnDiag2++;
		}
		if (i == 0) break;
	}

	prevOnDiag = false; counter = 0;
	for (i = n - 1; i >= 0; --i) {
		if (prevOnDiag == false and reverse_onDiag[i] == true) 
			diagStart = i;
		if (prevOnDiag == true and reverse_onDiag[i] == false) {
			if (diagStart - i + 1 < fopts.DiagMinCluster) { // [diagStart, i]
				for (j = i; j <= diagStart; ++j) { reverse_onDiag[j] = false; }
			}
			else { reverse_onDiag[i] = true; }
			counter++;
		}
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
			c++;
		}
	}
	hits.resize(c);
};

















