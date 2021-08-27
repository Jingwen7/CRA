#include "genome.h"
#include "index.h"
#include "option.h"
#include "input.h"
#include "rfpriv.h"
#include "hit.h"
#include "cluster.h"
#include "breakpoint.h"
#include "sample.h"
#include "graph.h"
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <errno.h>
#include <zlib.h>
#include <set>
#include <vector>


int main(int argc, char *argv[])
{
	if (argc < 2) {
		fprintf(stderr, "Usage: rf target.fa\n");
		return 1;
	}

	idxopt_t siopt(20, 1); // cannot be too small, 15 is a good start
	idxopt_t liopt(80, 1);
	fragopt_t fopts;

	// read input FASTA file
	Input reader;
	reader.Initialize(argv[1]);
	Genome genome;
	genome.read(argv[1]);

	// read index
	// if not available, generate index
	string sidxFile = string(argv[1]) + ".sidx";
	string lidxFile = string(argv[1]) + ".lidx";
	idx_t mi_s(siopt, genome);
	idx_t mi_l(liopt, genome);
	// TODO(Jingwen): seems the index is not readable
	if (mi_s.readIndex(sidxFile) == 1) {
		cerr << "k: " << mi_s.k << endl;
		mi_s.idx_gen(fopts);
		mi_s.idx_sort(fopts);
		mi_s.storeIndex(sidxFile);
	}
	if (mi_l.readIndex(lidxFile) == 1) {
		cerr << "k: " << mi_l.k << endl;
		mi_l.idx_gen(fopts);
		mi_l.idx_sort(fopts);
		mi_l.storeIndex(lidxFile);
	}
	cerr << "finish indexing!" << endl;

	// find small & dense hits for each sample
	uint32_t i, j;
	uint32_t n = genome.getSize();

	vector<sample> samples(genome.getSize());
	vector<vector<acrosample>> acrosamples(n);
	graph superGraph, hyperGraph;


	FindBasicRepeatUnit(samples, acrosamples, genome, mi_s, mi_l, siopt, liopt, fopts);

	// construct the graph
	// for each samples[i], find representative intervals as vertices
	// add edge between any pair of intervals that have matches
	// vector<uint32_t> clusterPivots;

	graphConstruction (superGraph, n, samples, acrosamples, fopts, 1);
	
	// project based on sparse_clusts, and cluster hyperbreakpoints based on initial breakponts
	FindsmallerRepeatUnit (samples, acrosamples, genome, mi_s, mi_l, siopt, liopt, fopts);
	graphConstruction (hyperGraph, n, samples, acrosamples, fopts, 0);

	return 0;
}
