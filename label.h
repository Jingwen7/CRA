#ifndef LABEL_H_
#define LABEL_H_

#include <vector>
#include <string>

class color_assign
{
	uint32_t colors;
	vector<uint32_t> c_vec;
};

typedef vector<color_assign> color_assigns;

class fg
{
public:
	static enum fgtype { INS, DEL, NONE};
	vector<uint32_t> len;
	vector<fgtype> type;
	vector<color_assigns> casn;
};

class sample
{
public:
	string name;
	string chrom;
	int cn;
	uint32_t clen;
	vector<fg> fgs; // store length/type/color_assignment of each fragment
	void unify (vector<copy_label> &cl, );
};

class 


#endif