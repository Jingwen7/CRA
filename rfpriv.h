#ifndef RFPRIV_H
#define RFPRIV_H

#include <assert.h>
#include <boost/multiprecision/cpp_int.hpp>
#include <map>
#include <vector>

using namespace boost::multiprecision;
using namespace std;

#define BITLIMIT 320
#define KMERBITS 256
#define RIDBITS 32
#define POSBITS 32
#define UINT256_MAX ~uint256_t(0)
#define RF_IDX_MAGIC "RFI\2"

// static for_mask = (((uint256_t)1 << (256 - 1)) - 1) << 1; // 111...1110
// static rev_mask = ~formask;

// emulate 288-bit integers and arrays
typedef struct { uint256_t x; uint64_t y; } idx_320_t; // ((<=254) minimizer + (1)direction: 0/1 + (64) rid + position)) 
typedef map<uint256_t, vector<uint32_t>> kmerhash; // (<=254) minimizer + (1)direction: 0/1, 32 pos
enum pointType { s, e };
typedef tuple<pointType, uint32_t, uint32_t> pointInfo; // point type, coord, cluster_index
typedef struct { uint32_t s; uint32_t e; uint32_t code; uint32_t cidx;} fragmatch; // start, end, code, cluster_index 
typedef pair<uint32_t, uint32_t> intv;
// void sketch(const char *str, int len, int w, int k, uint32_t rid, vector<idx_320_t> &idxv, uint32_t totalLen, bool is_hpc);
void sketch(char *seq, uint32_t seqLen, uint32_t k, uint32_t w, uint32_t rid, vector<idx_320_t> &idxv, bool debug);


#endif
