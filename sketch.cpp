#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include "rfpriv.h"

// static idx_256_t rev_mask = 1 << (256 - 1);
// static idx_256_t for_mask = ~ (rev_mask);	

unsigned char seq_nt4_table[256] = { // Table to change "ACGTN" to 01234
	0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};


// *
//  * Find asymmetric (w,k)-minimizers on a DNA sequence
//  *
//  * @param str    DNA sequence
//  * @param len    length of $str
//  * @param w      find a minimizer for every $w consecutive k-mers
//  * @param k      k-mer size
//  * @param rid    reference ID; will be copied to the output $p array
//  * @param idxv   minimizers
//  *               idxv[i].x = kMer<<1 | strand
//  *               idxv[i].y = rid<<32 | lastPos
//  *               where lastPos is the position of the last base of the i-th minimizer,
//  *               and strand indicates whether the minimizer comes from the top or the bottom strand.

 
void sketch(const char *str, int len, int w, int k, uint32_t rid, vector<idx_320_t> &idxv)
{
	// uint256_t shift1 = 2 * (k - 1);
	uint256_t kmer[2] = {(uint256_t) 0, (uint256_t) 0};
	// uint256_t formask = (1ULL << 2*k) - 1; // 
	int i, j, l, buf_pos, min_pos = 0;
	idx_320_t buf[256];
	idx_320_t min = { UINT256_MAX, UINT64_MAX };

	assert(len > 0 && (w > 0 && w < 256) && (k > 0 && k <= 127)); // 254 bits for k-mer; 
	memset(buf, 0xff, w * 16);
	idxv.resize(len/w);
	// kv_resize(mm128_t, km, *p, p->n + len/w);

	for (i = l = buf_pos = min_pos = 0; i < len; ++i) {
		int c = seq_nt4_table[(uint8_t)str[i]];
		idx_320_t info = { UINT256_MAX, UINT64_MAX};
		if (c < 4) { // not an ambiguous base
			int z;
  			kmer[0] = kmer[0]<<2 | c;           // forward k-mer
			kmer[1] = (kmer[1]>>2) | ((3^c) << (2 * (k - 1))); // reverse k-mer
			if (kmer[0] == kmer[1]) continue; // skip "symmetric k-mers" as we don't know it strand
			z = kmer[0] < kmer[1]? 0 : 1; // strand
			++l;
			if (l >= k) {
				info.x = kmer[z]<<1 | z;
				info.y = (uint64_t)rid<<32 | (uint32_t)i;
			}
		} 
		else {
			l = 0;
		}
		buf[buf_pos] = info; // need to do this here as appropriate buf_pos and buf[buf_pos] are needed below
		if (l == w + k - 1 && min.x != UINT256_MAX) { // special case for the first window - because identical k-mers are not stored yet
			for (j = buf_pos + 1; j < w; ++j)
				if (min.x == buf[j].x && buf[j].y != min.y) idxv.push_back(buf[j]);
			for (j = 0; j < buf_pos; ++j)
				if (min.x == buf[j].x && buf[j].y != min.y) idxv.push_back(buf[j]);
		}
		if (info.x <= min.x) { // a new minimum; then write the old min
			if (l >= w + k && min.x != UINT256_MAX) idxv.push_back(min); 
			min = info, min_pos = buf_pos;
		} 
		else if (buf_pos == min_pos) { // old min has moved outside the window
			if (l >= w + k - 1 && min.x != UINT256_MAX) idxv.push_back(min); 
			for (j = buf_pos + 1, min.x = UINT256_MAX; j < w; ++j) // the two loops are necessary when there are identical k-mers
				if (min.x >= buf[j].x) min = buf[j], min_pos = j; // >= is important s.t. min is always the closest k-mer
			for (j = 0; j <= buf_pos; ++j)
				if (min.x >= buf[j].x) min = buf[j], min_pos = j;
			if (l >= w + k - 1 && min.x != UINT256_MAX) { // write identical k-mers
				for (j = buf_pos + 1; j < w; ++j) // these two loops make sure the output is sorted
					if (min.x == buf[j].x && min.y != buf[j].y) idxv.push_back(buf[j]); 
				for (j = 0; j <= buf_pos; ++j)
					if (min.x == buf[j].x && min.y != buf[j].y) idxv.push_back(buf[j]); 
			}
		}
		if (++buf_pos == w) buf_pos = 0;
	}
	if (min.x != UINT256_MAX)
		idxv.push_back(min);
}


