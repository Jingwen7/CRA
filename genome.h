#ifndef GENOME_H_
#define GENOME_H_

#include <map>
#include <vector>
#include <assert.h>
#include <iostream>
#include <fstream>
#include <istream>
#include <ctype.h>
#include <string>
#include <zlib.h>
#include "htslib/kseq.h"
#include "htslib/hts.h"
#include "input.h"

using namespace std;

// declare the type of file handler and the read() function  
KSEQ_INIT(gzFile, gzread);

class Genome 
{
 public:
 	vector<string> names; 
	vector<char*> seqs;
	vector<uint32_t> lens;
	map<string, int> nameMap;

	Genome () {}

	~Genome() 
	{
		for (int i = 0; i < seqs.size(); i++) {
			delete[] seqs[i];
			seqs[i] = NULL;
		}
	}

	void add(const char* name, uint64_t l) 
	{
		names.push_back(string(name));
		lens.push_back(l);
	}

	void resize(uint32_t size)
	{
		names.resize(size);
		seqs.resize(size);
		lens.resize(size);
	}

	int getSize() 
	{
		return names.size();
	}

	string *getName (int index) 
	{
		return &names[index];
	}

	uint32_t getLen (int index) 
	{
		return lens[index];
	}

	char *getSeq (int index) 
	{
		return seqs[index];
	}

	char *subSeq(uint32_t index, const string &s) 
	{
		assert(nameMap[s] < getSize());
		return &seqs[nameMap[s]][index];
	}

	void read(const char* genomefile) 
	{
		gzFile fp = gzopen(genomefile, "r");
		// FILE* fp = fopen(genomefile.c_str(), "r");
		kseq_t *ks = kseq_init(fp); // initialize seq 
		int i = 0;
		while (kseq_read(ks) >= 0) { // each kseq_read() call reads one query sequence
			char *seq = new char[ks->seq.l];
			for (int j=0; j< ks->seq.l; j++) { seq[j] = toupper(ks->seq.s[j]); } // Convert lowercase letter to uppercase
			seqs.push_back(seq);
			lens.push_back(ks->seq.l);
			names.push_back(ks->name.s);
			nameMap[ks->name.s] = i;
			i++;
		}
		kseq_destroy(ks);
		gzclose(fp); // close the file handler
	}

	int getIndex(const string &s) 
	{
		if (nameMap.count(s) > 0) {
			cerr << s << " not in nameMap" << endl;
			exit(1);			
		}
		else { return nameMap[s]; }
	}

};

#endif
