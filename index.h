#ifndef INDEX_H_
#define INDEX_H_

#include <numeric>     /*iota */
#include <cmath>       /* ceil */
#include <utility>      // std::pair, std::make_pair
#include <map>
#include "rfpriv.h"
#include "genome.h"
#include "option.h"

class idx_t {
public: 
	uint32_t k, w;
	uint32_t nseq;
	int dnkmer;
	int sgnkmer;
	vector<idx_320_t> idxv; // ((<=254) minimizer + (1)direction: 0/1 + (64) rid + position)) vector
	vector<kmerhash *> seqh; // hash table (global pos)
	Genome *genome;	// seq, name, len
	idxopt_t *opt;

	idx_t (idxopt_t &iopt, Genome &g) : opt(&iopt), genome(&g) 
	{
		k = iopt.k;
		w = iopt.w;
		nseq = genome->getSize();
		seqh.resize(nseq);
		dnkmer = sgnkmer = 0;
	};

	~idx_t () {};

	void idx_init (const uint32_t K, const uint32_t W, const uint32_t N, const int d, const int s)
	{
		k = K;
		w = W;
		nseq = N;
		dnkmer = d;
		sgnkmer = s;
	}

	void idx_gen ();

	void idx_sort ();
	
	void idx_stats ();

	int readIndex (string fn);

	int storeIndex (string fn);
};

class sortPreBitsOp 
{
public:
    static uint256_t pmask;
    sortPreBitsOp (int bits) 
    {
    	pmask = ((uint256_t)1 << (bits - 1) - 1) << 1; // 111...1110
    }
    ~sortPreBitsOp () {};
    bool operator() (const idx_320_t &a, const idx_320_t &b) 
    {
    	return (a.x & pmask) < (b.x & pmask);
    }
};

class sortSufBitsOp 
{
public:
    static uint64_t smask;
    sortSufBitsOp (int bits) 
    {
    	smask = (uint64_t)1 << bits - 1; // 0000....1111
    }
    ~sortSufBitsOp () {};
    bool operator() (const idx_320_t &a, const idx_320_t &b)  
    {
		return (a.y & smask) < (b.y & smask);
    }
}; 
#endif
