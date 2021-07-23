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

// emulate 288-bit integers and arrays
typedef struct { uint256_t x; uint64_t y; } idx_320_t;
typedef map<uint256_t, vector<uint32_t>> kmerhash; // (<=254) minimizer + (1)direction: 0/1, 32 pos
typedef struct { uint32_t x; uint32_t y; } hit;

void sketch(const char *str, int len, int w, int k, uint32_t rid, vector<idx_320_t> &idxv);

void hit(const string &s_a, const string &s_b, const idx_t &mi_s);

#endif
