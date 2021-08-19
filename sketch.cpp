#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include "rfpriv.h"
#include "index.h"
	
uint256_t for_mask = (((uint256_t)1 << (256 - 1)) - 1) << 1; // 111...1110
uint256_t rev_mask = ~for_mask; // 0000....0001;

// debug
void checkIfKmerCorrect (idx_320_t minimizer, uint32_t k, char *seq) {
	// debug
	bool strand = (bool) (minimizer.x & rev_mask);
	minimizer.x = (minimizer.x & for_mask) >> 1;
	checkMinimizerMatch(minimizer.x, k, seq, minimizer.y, strand);
}

void printKmer (idx_320_t minimizer, uint32_t k, char *seq) {
	// debug
	bool strand = (bool) (minimizer.x & rev_mask);
	minimizer.x = ((minimizer.x & for_mask) >> 1);
	printMinimizerMatch(minimizer.x, k, seq, minimizer.y, strand);
}

unsigned int seq_nt4_table[256] = { // Table to change "ACGTN" to 01234
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

char BinaryTont4_table[4] = { // Table to change 01234 "ACGTN" to
	'A', 'C', 'G', 'T'
};


void InitMask(idx_320_t &m, uint32_t k) {
	// m.x = 001111111
	m.x = (uint256_t) 0;
	for (uint32_t i = 0; i < k; i++) {
		m.x <<= 2;
		m.x += 3;
	}
}

void StoreKmer(char *seq, uint32_t pos, uint32_t k, idx_320_t &t) {
	t.x = (uint256_t) 0;
	for (uint32_t p = pos; p <= pos + k - 1 ; p++) {
		t.x <<= 2;
		t.x += seq_nt4_table[(uint8_t)seq[p]];
	}
	t.x <<= 1; // shift one left for last strand bit
}

void ShiftOne(char *seq, uint32_t pos, const idx_320_t &mask, idx_320_t &t) {
	t.x >>= 1; // shift one for last strand bit
	t.x = (t.x << 2) & mask.x;
	t.x += seq_nt4_table[(uint8_t)seq[pos]];	
	t.x <<= 1; // shift one for last strand bit
}

void ShiftOneRC(char *seq, uint32_t pos, uint32_t k, idx_320_t &t) {
	t.x >>= 1; // shift one for last strand bit	
	int c = seq_nt4_table[(uint8_t)seq[pos]];	
	t.x >>= 2;
	t.x |= ((3 ^ (uint256_t) c) << (2 * (k - 1)));
	t.x <<= 1; // shift one for last strand bit
	t.x |= rev_mask;
}

void KmerRC(idx_320_t a, idx_320_t &b, uint32_t k, char * seq) {
	a.x >>= 1; // shift one for last strand bit
	uint256_t nucMask = 3;
	uint256_t least;
	b.x = (uint256_t) 0;
	for (uint32_t i = 0; i < k; i++) {
		least = (a.x & nucMask) ^ nucMask;
		a.x >>= 2;
		b.x <<= 2;
		b.x += least;
	}
	b.x <<= 1; // shift one for last strand bit
	b.x |= rev_mask;
}

void sketch(char *seq, uint32_t seqLen, uint32_t k, uint32_t w, uint32_t rid, vector<idx_320_t> &idxv, bool debug) {
	//
	// Initialize first.
	//
	if (seqLen < k) {
		return;
	}
	idx_320_t cur, curRC, can;
	bool strand;
	uint32_t minPos = 0, p = 0;
	int windowSpan = w + k - 1;
	idx_320_t m;
	//
	// Skip N's as a start
	//
	InitMask(m, k); 
	int nextValidWindowEnd = 0;
	int nextValidWindowStart = 0;
	bool valid = false;
	if (seqLen < windowSpan) return; 
	while (nextValidWindowStart < seqLen - windowSpan and !valid) {
		valid = true;
		for (int n = nextValidWindowStart; valid and n < nextValidWindowStart + windowSpan; n++) {
			if (seqLen < n) return;
			if (seq_nt4_table[(uint8_t)seq[n]] > 3) {
				nextValidWindowStart = n + 1;
				valid = false;
			}
		}
	}
	// all n
	if (valid == false) {
		return;
	}
	nextValidWindowEnd = nextValidWindowStart + windowSpan;

	StoreKmer(seq, p, k, cur);
	KmerRC(cur, curRC, k, seq);

	//
	// Initialize the first minimzer.
	//
	if ((cur.x & for_mask) < (curRC.x & for_mask)) can.x = cur.x;
	else can.x = curRC.x; 

	idx_320_t activeKmer, curKmer, insertKmer;
	activeKmer.x = can.x;
	activeKmer.y = 0;
	vector<idx_320_t> curKmers(w);
	curKmers[0] = activeKmer;

	// 
	// Find the active minimizer in this window
	//
	int nKmer = 1;
	int nSkipped = 0;

	for (p = 1; p < w && p < seqLen - k + 1 ; p++) {

		ShiftOne(seq, p + k - 1, m, cur);
		ShiftOneRC(seq, p + k - 1, k, curRC);
		/*
		Tuple test, testrc;
		StoreKmer(seq->seq.s, p, k, test);
		KmerRC(test, testrc, k);
		assert(test == cur);
		assert(testrc == curRC);
		*/
		curKmer.y = p;
		if ((cur.x & for_mask) < (curRC.x & for_mask)) curKmer.x = cur.x;
		else curKmer.x = curRC.x;
		// if (debug) checkIfKmerCorrect(curKmer, k, seq);

		if (curKmer.x < activeKmer.x) {  
			activeKmer.x = curKmer.x;
			activeKmer.y = p;
		}	
		curKmers[p % w] = curKmer;
	}
	//
	// Only store the first minimizer if the first window starts at the beginning of the sequence.
	//
	if (nextValidWindowEnd == windowSpan) {
		insertKmer = activeKmer;
		// debug
		uint256_t minimizer = (insertKmer.x & for_mask) >> 1;
		strand = (bool) (insertKmer.x & rev_mask);
		checkMinimizerMatch(minimizer, k, seq, insertKmer.y, strand);					

		insertKmer.y += ((uint64_t)rid << 32);
		idxv.push_back(insertKmer);
	}
	// Now scan the chromosome
	for (p = w; p < seqLen - k + 1; p++) {
		// If the next valid window ends at the next nucleotide, check to see if 
		// it is a valid window (no N's). If so, bump by one.
		// Otherwise, search for the next valid window end.
		if (nextValidWindowEnd == p + k - 1)  {
			if ((uint8_t)seq_nt4_table[seq[p+k-1]] <= 3) {
				nextValidWindowEnd++;
			}
			else {
				nextValidWindowStart = p + k;
				valid = false;		
				if (seqLen < windowSpan) return; 	
				while (nextValidWindowStart < seqLen - windowSpan and not valid) {
					valid = true;
					for (int n = nextValidWindowStart; valid and n < nextValidWindowStart+windowSpan; n++) {
						if ((uint8_t)seq_nt4_table[seq[n]] > 3) {
							nextValidWindowStart = n + 1;
							valid = false;
						}
					}			
				}
				// all n
				if (valid == false) {
					return;						
				}
				nextValidWindowEnd = nextValidWindowStart + windowSpan;
			}
		}		

		ShiftOne(seq, p+k-1, m, cur);
		ShiftOneRC(seq, p+k-1, k, curRC);

// #ifdef _TESTING_
// 		uint256_t test, testrc;
// 		StoreKmer(seq, p, k, test);
// 		KmerRC(test, testrc, k);

// 		assert(test.x == cur.x);
// 		assert(testrc.x == curRC.x);
// #endif
		if ((cur.x & for_mask) < (curRC.x & for_mask)) curKmer.x = cur.x;
		else curKmer.x = curRC.x; 
		curKmer.y = p;
		// if (debug) checkIfKmerCorrect(curKmer, k, seq);
		curKmers[p % w] = curKmer;
		if (p - w >= activeKmer.y) {
			activeKmer = curKmers[0];
			for (uint32_t j = 1; j < w; j++) {
				if ((curKmers[j].x & for_mask) < (activeKmer.x & for_mask)) { 
					activeKmer = curKmers[j];
				}		
			}
			if (nextValidWindowEnd == p + k) {
				insertKmer = activeKmer;
				// debug
				// if (debug) {
				// 	uint256_t minimizer = (insertKmer.x & for_mask) >> 1;
				// 	strand = (bool) (insertKmer.x & rev_mask);
				// 	checkMinimizerMatch(minimizer, k, seq, insertKmer.y, strand);						
				// }

				insertKmer.y += ((uint64_t)rid << 32);
				idxv.push_back(insertKmer);			
				nKmer+=1;
			}
			else {
				nSkipped++;
			}
		}
		else {
			if ((curKmer.x & for_mask) < (activeKmer.x & for_mask)) { //TODO(Jingwen)
				activeKmer = curKmer;
				if (nextValidWindowEnd == p+k) {
					insertKmer = activeKmer;
					// // debug
					// if (debug) {
					// 	uint256_t minimizer = (insertKmer.x & for_mask) >> 1;
					// 	strand = (bool) (insertKmer.x & rev_mask);
					// 	checkMinimizerMatch(minimizer, k, seq, insertKmer.y, strand);						
					// }	

					insertKmer.y += ((uint64_t)rid << 32);
					idxv.push_back(insertKmer);	
					nKmer++;
				}
				else {
					nSkipped++;
				}
			}		
		}		
		if (p + 1 % 10000 == 0) {
			cerr << p + 1 << "\t" << nKmer << "\t" << nSkipped << endl;
		}
	}

	for (uint32_t i = 0; i < idxv.size(); ++i) {
		// if (debug) checkIfKmerCorrect(curKmer, k, seq);
	}
}













































































// // typedef struct { // a simplified version of kdq
// // 	int front, count;
// // 	int a[32];
// // } tiny_queue_t;

// // static inline void tq_push(tiny_queue_t *q, int x)
// // {
// // 	q->a[((q->count++) + q->front) & 0x1f] = x;
// // }

// // static inline int tq_shift(tiny_queue_t *q)
// // {
// // 	int x;
// // 	if (q->count == 0) return -1;
// // 	x = q->a[q->front++];
// // 	q->front &= 0x1f;
// // 	--q->count;
// // 	return x;
// // }


// // *
// //  * Find asymmetric (w,k)-idxv on a DNA sequence
// //  *
// //  * @param str    DNA sequence
// //  * @param len    length of $str
// //  * @param w      find a minimizer for every $w consecutive k-mers
// //  * @param k      k-mer size
// //  * @param rid    reference ID; will be copied to the output $p array
// //  * @param idxv   idxv
// //  *               idxv[i].x = kMer<<1 | strand
// //  *               idxv[i].y = rid<<32 | lastPos
// //  *               where lastPos is the position of the last base of the i-th minimizer,
// //  *               and strand indicates whether the minimizer comes from the top or the bottom strand.

 
// void sketch(const char *str, int len, int w, int k, uint32_t rid, vector<idx_320_t> &idxv, uint32_t totalLen, bool is_hpc)
// {
// 	// uint256_t shift1 = 2 * (k - 1);
// 	uint256_t kmer[2] = {(uint256_t) 0, (uint256_t) 0};
// 	// uint256_t formask = (1ULL << 2*k) - 1; // 
// 	int i, j, l, buf_pos, min_pos = 0, kmer_span = 0;
// 	idx_320_t buf[256];
// 	idx_320_t min = { UINT256_MAX, UINT64_MAX };
// 	tiny_queue_t tq;

// 	assert(len > 0 && (w > 0 && w < 256) && (k > 0 && k <= 127)); // 254 bits for k-mer; 
// 	memset(buf, 0xff, w * 16);
// 	memset(&tq, 0, sizeof(tiny_queue_t));
// 	idxv.resize(totalLen/w);

// 	for (i = l = buf_pos = min_pos = 0; i < len; ++i) {
// 		int c = seq_nt4_table[(uint8_t)str[i]];
// 		idx_320_t info = { UINT256_MAX, UINT64_MAX};
// 		if (c < 4) { // not an ambiguous base
// 			int z;
// 			if (is_hpc) {
// 				int skip_len = 1;
// 				if (i + 1 < len && seq_nt4_table[(uint8_t)str[i + 1]] == c) {
// 					for (skip_len = 2; i + skip_len < len; ++skip_len)
// 						if (seq_nt4_table[(uint8_t)str[i + skip_len]] != c)
// 							break;
// 					i += skip_len - 1; // put $i at the end of the current homopolymer run
// 				}
// 				tq_push(&tq, skip_len);
// 				kmer_span += skip_len;
// 				if (tq.count > k) kmer_span -= tq_shift(&tq);
// 			} 
// 			else {
// 				kmer_span = l + 1 < k? l + 1 : k;
// 			}
//   			kmer[0] = kmer[0]<<2 | c;           // forward k-mer
// 			kmer[1] = (kmer[1]>>2) | ((3^c) << (2 * (k - 1))); // reverse k-mer
// 			if (kmer[0] == kmer[1]) continue; // skip "symmetric k-mers" as we don't know it strand
// 			z = kmer[0] < kmer[1]? 0 : 1; // strand
// 			++l;
// 			if (l >= k and kmer_span < 256) {
// 				info.x = kmer[z]<<1 | z;
// 				info.y = (uint64_t)rid<<32 | (uint32_t)i;
// 			}
// 		} 
// 		else {
// 			l = 0, tq.count = tq.front = 0, kmer_span = 0;
// 		}
// 		buf[buf_pos] = info; // need to do this here as appropriate buf_pos and buf[buf_pos] are needed below
// 		if (l == w + k - 1 && min.x != UINT256_MAX) { // special case for the first window - because identical k-mers are not stored yet
// 			for (j = buf_pos + 1; j < w; ++j)
// 				if (min.x == buf[j].x && buf[j].y != min.y) idxv.push_back(buf[j]);
// 			for (j = 0; j < buf_pos; ++j)
// 				if (min.x == buf[j].x && buf[j].y != min.y) idxv.push_back(buf[j]);
// 		}
// 		if (info.x <= min.x) { // a new minimum; then write the old min
// 			if (l >= w + k && min.x != UINT256_MAX) idxv.push_back(min); 
// 			min = info, min_pos = buf_pos;
// 		} 
// 		else if (buf_pos == min_pos) { // old min has moved outside the window
// 			if (l >= w + k - 1 && min.x != UINT256_MAX) idxv.push_back(min); 
// 			for (j = buf_pos + 1, min.x = UINT256_MAX; j < w; ++j) // the two loops are necessary when there are identical k-mers
// 				if (min.x >= buf[j].x) min = buf[j], min_pos = j; // >= is important s.x. min is always the closest k-mer
// 			for (j = 0; j <= buf_pos; ++j)
// 				if (min.x >= buf[j].x) min = buf[j], min_pos = j;
// 			if (l >= w + k - 1 && min.x != UINT256_MAX) { // write identical k-mers
// 				for (j = buf_pos + 1; j < w; ++j) // these two loops make sure the output is sorted
// 					if (min.x == buf[j].x && min.y != buf[j].y) idxv.push_back(buf[j]); 
// 				for (j = 0; j <= buf_pos; ++j)
// 					if (min.x == buf[j].x && min.y != buf[j].y) idxv.push_back(buf[j]); 
// 			}
// 		}
// 		if (++buf_pos == w) buf_pos = 0;
// 	}
// 	if (min.x != UINT256_MAX)
// 		idxv.push_back(min);
// }


