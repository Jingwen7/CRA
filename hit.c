#include rfpriv.h
#include <vector>
#include <map>

inline void find_match (vector<hit> &hits, const uint256_t &mi, const vector<uint32_t> &pos, const kmerhash &khash, bool self)
{
	if (self) {
	// (TODO) Jingwen: add code to only maintain half of the hits
	}
	else {
		auto ot;
		ot = khash.find(mi);
		if (ot != khash.end) {
			uint32_t i, j;
			for (i = 0; i < ot->second.size(); ++i) {
				for (j = 0; j < pos.size(); ++j) {
					hit temp;
					temp.x = pos[j]; temp.y = (ot->second)[i];
					hits.push_back(temp);				
				}
			}
		}		
	}
}

void hit(const uint32_t a, const uint32_t b, const idx_t &mi, vector<hit> &fhits, vector<hit> &rhits) 
{
	bool self = 1;
	if (a != b) self = 0; 

	for (auto it = mi.seqh[a]->begin(); it! = mi.seqh[a]->end(); ++it) {
		it->first // minimizer
		it->second // vector<pos>

		uint256_t rmask, fmask, z, f, r;
		rmask = (uint256_t) 1; // 00000...1
		fmask = (((uint256_t)1 << 255) - 1) << 1; // 111...1110
		z = ((rmask & it->first) == 1) ? 1 : 0; // direction
		f = (it->first & fmask);
		r = (it->first | rmask);

		if (z ^ 0 == 1) // reverse
			find_match(rhits, f, mi.seqh[a]->second, *mi.seqh[b]);
		else if (z ^ 1 == 1) // reverse
			find_match(rhits, r, mi.seqh[a]->second, *mi.seqh[b]);
		else if (z ^ 0 == 0) // forward
			find_match(fhits, f, mi.seqh[a]->second, *mi.seqh[b]);
		else if (z ^ 1 == 0) // forward
			find_match(fhits, r, mi.seqh[a]->second, *mi.seqh[b]);
	}
}
