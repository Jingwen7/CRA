#include <map>
#include <fstream>
#include <iostream>
#include "index.h"
#include "rfpriv.h"

using namespace std;

uint256_t sortPreBitsOp::pmask;
uint64_t sortSufBitsOp::smask;


void idx_t::idx_gen () 
{
	for (int i = 0; i < genome->getSize(); i++) {
		sketch(genome->getSeq(i), genome->getLen(i), w, k, i, idxv);
	}
};

void idx_t::idx_sort () 
{
	// sort by 256 bits (idx_320_t.x)
	sortPreBitsOp presorter(KMERBITS);
	sort(idxv.begin(), idxv.end(), presorter);

	// create hash table
	int s = 0;
	uint64_t pmask = ((uint64_t)1 << POSBITS - 1); // 0000...1111
	uint64_t rmask = ~pmask; // 1111...0000
	dnkmer = 0; sgnkmer = 0;

	for (int i = 1; i <= idxv.size(); i++) {
		if (i < idxv.size() and (idxv[i].x & presorter.pmask) == (idxv[i - 1].x & presorter.pmask)) continue;
		else {
			dnkmer++;
			if (i - s > 1) {
				// sort by pos
				sortSufBitsOp sufsorter(POSBITS);
				if (i < idxv.size()) sort(idxv.begin() + s, idxv.begin() + i, sufsorter);
				else sort(idxv.begin() + s, idxv.end(), sufsorter);
			}
			else sgnkmer++;
			for (int j = s; j < i; j++) {
				uint32_t rid = (uint32_t) ((idxv[j].y & rmask) >> POSBITS);
				uint32_t pos = (uint32_t) (idxv[j].y & pmask);
				// string name = genome->getName(rid);
				// f = idxv.x & presorter.mask; // forward strand
				// r = idxv.x & presorter.mask + 1; // backward strand
				// z = idxv.x == f? 0 : 1; // direction
				// seqh[rid] 
				if (seqh[rid]->empty()) (*seqh[rid])[idxv[j].x].push_back(pos);
				else {
					kmerhash::iterator t = seqh[rid]->find(idxv[j].x);
					if (t != seqh[rid]->end()) (*seqh[rid])[idxv[j].x].push_back(pos);
				}
			}
			s = i;
		}
	}
	idxv.clear();
}

void idx_t::idx_stats () 
{
	uint32_t i;
	uint64_t sum = 0, len = 0;
	fprintf(stderr, "[M::%s] kmer size: %d; skip: %d; #seq: %d\n", __func__, k, w, nseq);
	for (i = 0; i < nseq; ++i) {
		len += genome->getLen(i);
	}
	for (i = 0; i < nseq; ++i) {
		if (!seqh[i]->empty()) {
			for (kmerhash::iterator t = seqh[i]->begin(); t != seqh[i]->end(); ++t) {
				sum += (t->second).size();
			}
		}
	}
	fprintf(stderr, "[M::%s] distinct minimizers: %d (%.2f%% are singletons); average occurrences: %.3lf; average spacing: %.3lf; total length: %ld\n",
			__func__, dnkmer, 100.0 * sgnkmer / dnkmer, (double)sum / dnkmer, (double)len / sum, (long)len);	
}

int idx_t::storeIndex (string fn) {
	ofstream fout(fn.c_str(), ios::out|ios::binary);
	fout.write(RF_IDX_MAGIC, 4);
	uint32_t x[3];
	x[0] = k, x[1] = w, x[2] = nseq;
	fout.write((char*) x, sizeof(uint32_t) * 3);
	fout.write((char*) &dnkmer, sizeof(int));
	fout.write((char*) &sgnkmer, sizeof(int));

	uint32_t i;
	uint32_t nkmer;
	for (i = 0; i < nseq; ++i) {
		// write seq_name, seq_len
		size_t size = (genome->getName(i))->size();
		fout.write((char *) &size, sizeof(size));
		fout.write((genome->getName(i))->c_str(), size);
		uint32_t len = genome->getLen(i);
		fout.write((char*) &len, sizeof(uint32_t));

		// write hash table
		if (!seqh[i]->empty()) {
			nkmer = seqh[i]->size();
			fout.write((char*) &nkmer, sizeof(uint32_t));
			for (kmerhash::iterator t = seqh[i]->begin(); t != seqh[i]->end(); ++t) {
				fout.write((char*) &(t->first), sizeof(uint256_t));
				size = (t->second).size();
				fout.write((char*) &size, sizeof(size));
				fout.write((char*) &(t->second), sizeof(uint32_t)*(t->second).size());
			}				
		}
		else {
			// uint32_t 
			// fout.write();
		}
	}

	fout.close();
	return 0;
}

// magic word, k, w, nseq, dnkmer, sgnkmer,(seqname, seqlen, (nkmer, (kmer, pos_vec))
int idx_t::readIndex (string fn) 
{
	ifstream fin(fn.c_str(), ios::in|ios::binary);
	if (fin.good() == false or fin.eof()) {
    	// cerr << "Cannot open target " << filename << endl;
    	// exit(1);
		return 1;
	}
	char magic[3];
	uint32_t x[3];
	if (!fin.read(magic, 4)) return 1;
	if (strncmp(magic, RF_IDX_MAGIC, 4) != 0) return 1;
	if (!fin.read((char*)x, sizeof(uint32_t) * 3)) return 1;

	int d, s;
	if (!fin.read((char*) &d, sizeof(int))) return 1;
	if (!fin.read((char*) &s, sizeof(int))) return 1;
	idx_init(x[0], x[1], x[2], d, s);
	genome->resize(nseq);

	uint32_t i;
	for (i = 0; i < nseq; ++i) {
		size_t len;
		string name;
		if (!fin.read((char *)&len, sizeof(size_t))) return 1;
		char* temp = new char[len + 1];
		if (!fin.read(temp, len)) return 1;
		temp[len] = '\0';
		name = temp;
		delete [] temp;

		uint32_t seqlen;
		if (!fin.read((char*) &seqlen, sizeof(uint32_t))) return 1;

		(genome->names)[i] = name;
		(genome->lens)[i] = seqlen;
		(genome->nameMap)[name] = i;

		// read kmer and pos vector
		uint32_t nkmer;
		if (!fin.read((char*) &nkmer, sizeof(uint32_t))) return 1;
		uint32_t j;
		for (j = 0; j < nkmer; ++j) {
			uint256_t m;
			vector<uint32_t> pos;
			if (!fin.read((char*) &m, sizeof(uint256_t))) return 1;
			if (!fin.read((char*) &len, sizeof(size_t))) return 1;
			pos.resize(len);
			if (!fin.read((char*) &pos[0], sizeof(uint32_t) * len)); return 1;
			(*seqh[i])[m] = pos;
		}
	}
	return 0;
}