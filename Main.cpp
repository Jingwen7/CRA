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

	idxopt_t siopt(20, 1);
	idxopt_t liopt(50, 1);
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

	// Get the breakpoints for each sample dense_clusts
	vector<sample> samples(genome.getSize());
	for (i = 0; i < n; ++i) {
		cerr << "process sample: " << *genome.getName(i) << " " << i << endl;
		samples[i].init(i, genome.getName(i));
		samples[i].process(samples, i, i, mi_s, mi_l, fopts, siopt, liopt, 1, 1);
		samples[i].dump(genome.getName(i), fopts);
	}

	// get the dense_clusts for across-sample 
	vector<vector<acrosample>> acrosamples(n);
	for (i = 0; i < n; ++i) {
		acrosamples[i].resize(n - i);
		for (j = i + 1; j < n; ++j) {
			acrosamples[i][j - i - 1].init(i, j, &samples[i], &samples[j]);
			acrosamples[i][j - i - 1].across_process(samples, i, j, mi_s, mi_l, fopts, siopt, liopt, 0, 1);
			acrosamples[i][j - i - 1].decideStrand();
			// acrosamples[i][j - i - 1].trimclusters(1);
			if (fopts.debug)
				acrosamples[i][j - i - 1].dump(genome.getName(i), genome.getName(j), fopts, genome.getLen(j));
		}
	}	

	// project every sample to samples[i]
	if (n > 1) {
		for (i = 0; i < n; ++i) {
			for (j = 0; j < n; ++j) {
				// project samples[j] to samples[i]
				if (i == j) continue;
				if (j > i) { // samples[j] is y-axis
					project(samples[j], samples[i], acrosamples[i][j - i - 1].dense_clusts, acrosamples[i][j - i - 1].sparse_clusts, fopts, 0, 1);
				}
				else { // samples[j] is x-axis
					project(samples[j], samples[i], acrosamples[j][i - j - 1].dense_clusts, acrosamples[j][i - j - 1].sparse_clusts, fopts, 1, 1);
				}
			}
		}		
	}

	vector<uint32_t> clusterPivots;

	// for (i = 0; i < n; ++i) {
	// 	clusterPivots.clear();
	// 	samples[i].clusterBreakpoints(fopts, clusterPivots);
	// }

	// project every sample to samples[i]
	if (n > 1) {
		for (i = 0; i < n; ++i) {
			for (j = 0; j < n; ++j) {
				// project samples[j] to samples[i]
				if (i == j) continue;
				if (j < i) { // samples[i] is y-axis
					project(samples[i], samples[j], acrosamples[j][i - j - 1].dense_clusts, acrosamples[j][i - j - 1].sparse_clusts, fopts, 0, 1);
				}
				else { // samples[i] is x-axis
					project(samples[i], samples[j], acrosamples[i][j - i - 1].dense_clusts, acrosamples[i][j - i - 1].sparse_clusts, fopts, 1, 1);
				}
			}
		}		
	}

	for (i = 0; i < n; ++i) {
		for (j = i + 1; j < n; ++j) {
			if (fopts.debug)
				acrosamples[i][j - i - 1].dump(genome.getName(i), genome.getName(j), fopts, genome.getLen(j));
		}
	}	
	// for (i = 0; i < n; ++i) {
	// 	clusterPivots.clear();
	// 	samples[i].clusterBreakpoints(fopts, clusterPivots);
	// }

	// project self again
	for (i = 0; i < n; ++i) {
		selfproject (samples[i], samples[i].dense_clusts, samples[i].sparse_clusts, fopts, 1);
		clusterPivots.clear();
		samples[i].clusterBreakpoints(fopts, clusterPivots);
		selfproject (samples[i], samples[i].dense_clusts, samples[i].sparse_clusts, fopts, 1);
	}
	// if (fopts.debug) {
	// 	for (i = 0; i < n; ++i) {
	// 		samples[i].dump(genome.getName(i), fopts);
	// 		for (j = i + 1; j < n; ++j) {
	// 			acrosamples[i][j - i - 1].dump(genome.getName(i), genome.getName(j), fopts, genome.getLen(j));
	// 		}
	// 	}
	// }
	cerr << "finish projecting breakpoints" << endl;

	// construct the graph
	// for each samples[i], find representative intervals as vertices
	// add edge between any pair of intervals that have matches
	graph superGraph;

	// // decide relative strand for samples
	// samples[0].relative_strand = 0;
	// for (i = 1; i < n; ++i) {
	// 	if (acrosamples[i][j - 1].strand == 0) 
	// 		samples[i].relative_strand = 0;
	// 	else
	// 		samples[i].relative_strand = 1;
	// }
	for (i = 0; i < n; ++i) {
		// substract self-self dotplot clusters
		samples[i].substractClusters();
		clusterPivots.clear();
		samples[i].clusterBreakpoints(fopts, clusterPivots);

		// insert intervals to graph
		superGraph.insertInvt(clusterPivots, samples[i].idx, samples[i].s, samples[i].e);
	}

	superGraph.init();
	for (i = 0; i < n; ++i) 
		samples[i].unifyIntv(fopts, superGraph, 1);

	for (i = 0; i < n; ++i) {
		for (j = i + 1; j < n; ++j) {
			acrosamples[i][j - i - 1].unifyIntv(fopts, superGraph, 1);
		}
	}
	cerr << "finish constructing the graph" << endl;
	// find connected components in the superGraph
	superGraph.findConnetedComponents(1);
	dumpGraph(samples, superGraph, fopts);
	return 0;
}
