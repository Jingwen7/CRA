#include <map>
#include <fstream>
#include <iostream>
#include "index.h"
#include "rfpriv.h"

using namespace std;

uint256_t sortPreBitsOp::pmask;
uint64_t sortSufBitsOp::smask;

char binaryTont4_table[4] = { // Table to change 01234 "ACGTN" to
	'A', 'C', 'G', 'T'
};

void printK (uint256_t mi, uint32_t k) {
	char test[k];
	char * chararr = test;
	uint256_t mask = 3;
	uint32_t i;
	for (i = 0; i < k; ++i) {
		chararr[i] = binaryTont4_table[(int)(mi & mask)];
		cerr << chararr[i] ;
		mi >>= 2;
	}
	cerr << endl;	
}

void binaryTochar (uint256_t mi, char *chararr, uint32_t k, bool strand)
{
	uint256_t mask = 3;
	uint32_t i;
	for (i = 0; i < k; ++i) {
		// cerr << mi << " " << mask << " " << (mi & mask) << " " << (int)(mi & mask) << " "; 
		if (strand == 0) chararr[i] = binaryTont4_table[(int)(mi & mask)];
		else chararr[i] = binaryTont4_table[(int)((mi & mask) ^ mask)];
		mi >>= 2;
	}
	// cerr << endl;
}

void checkMinimizerMatch (uint256_t minimizer, uint32_t k, char * seq, uint32_t start, bool strand)
{
	char test[k];
	char * chararr = test;
	uint32_t d, l;
	binaryTochar(minimizer, chararr, k, strand);
	for (l = 0; l < k; ++l) {
		if (!strand) assert(chararr[k - 1 - l] == seq[start + l]);
		else assert(chararr[l] == seq[start + l]);
	}	
	// cerr << "minimizer: " << chararr << endl;	
}

void printMinimizerMatch (uint256_t minimizer, uint32_t k, char * seq, uint32_t start, bool strand)
{
	char test[k];
	char * chararr = test;
	uint32_t d, l;
	binaryTochar(minimizer, chararr, k, strand);
	for (l = 0; l < k; ++l) {
		cerr << chararr[l];
	}	
	cerr << " " << start << endl;	
}

void idx_t::idx_totalLen ()
{
	totalLen = accumulate(genome->lens.begin(), genome->lens.end(), 0);	
}

void idx_t::idx_gen (const fragopt_t &fopts) 
{
	for (int i = 0; i < genome->getSize(); i++) {
		sketch(genome->getSeq(i), genome->getLen(i), k, w, i, idxv, fopts.debug);
	}
	cerr << "idxv.size(): " << idxv.size() << endl;
};

void idx_t::idx_sort (const fragopt_t &fopts) 
{
	// sort by 256 bits (idx_320_t.x)
	sortPreBitsOp presorter(KMERBITS);
	sort(idxv.begin(), idxv.end(), presorter);
	// cerr << "presorter.pmask: " << presorter.pmask << endl;	
	ofstream fclust("minimizers.bed");

	//create hash table
	int s = 0;
	uint64_t pmask = (((uint64_t)1 << POSBITS) - 1); // 0000...1111
	uint64_t rmask = ~pmask; // 1111...0000
	dnkmer = 0; sgnkmer = 0;
	// cerr << "pmask: " << pmask << endl;
	// cerr << "rmask: " << rmask << endl;
	uint32_t i;
	for (i = 1; i <= idxv.size(); i++) {
		if (i < idxv.size() and (idxv[i].x & presorter.pmask) == (idxv[i - 1].x & presorter.pmask)) continue;
		else {
			dnkmer++;
			if (i - s > 1) {
				// sort by pos
				sortSufBitsOp sufsorter(POSBITS);
				if (i < idxv.size()) sort(idxv.begin() + s, idxv.begin() + i, sufsorter);
				else sort(idxv.begin() + s, idxv.end(), sufsorter);

				// // debug
				// if (fopts.debug) {
				// 	for (uint32_t m = s; m < i; ++m) {
				// 		fclust << idxv[m].x << "\t" << idxv[m].y << endl;
				// 	}					
				// }
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
				if (seqh[rid]->empty()) {
					seqh[rid]->insert(make_pair(idxv[j].x, vector<uint32_t>(1, pos))); //(*seqh[rid])[idxv[j].x].push_back(pos);
				}
				else {
					kmerhash::iterator t = seqh[rid]->find(idxv[j].x);
					if (t != seqh[rid]->end()) (*seqh[rid])[idxv[j].x].push_back(pos);
					else seqh[rid]->insert(make_pair(idxv[j].x, vector<uint32_t>(1, pos)));
				}

				//debug
				if (fopts.debug) {
					bool strand = (bool)(idxv[j].x & (~presorter.pmask));
					uint256_t minimizer = (idxv[j].x >> 1);
					checkMinimizerMatch(minimizer, k, genome->seqs[rid], pos, strand);						
				}
				
			}
			s = i;
		}
	}
	fclust.close();	
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
	char magic[4];
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

