#ifndef ALIGN_H_
#define ALIGN_H_

#include <vector>
#include <string>
#include 

class intv 
{
public:
	uint32_t a_s, a_e;
	uint32_t b_s, b_e;
};

class single_align_info
{
public:
	uint32_t idx_a, idx_b;
	vector<intv> fg;
};

typedef vector<single_align_info> single_align_infos;

class align_info 
{
public:
	uint32_t arid, brid;
	int na, nb; 
	vector<single_align_infos> ainfo;

	align_info (uint32_t Arid, uint32_t Brid, int Na, int Nb) : arid(Arid), brid(Brid), na(Na), nb(Nb)
	{
		ainfo.resize(na);
		for (int i = 0; i < na; ++i) { ainfo[i].resize(nb - i);}
	}
};

#endif