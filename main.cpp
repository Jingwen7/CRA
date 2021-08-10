#include "genome.h"
#include "index.h"
#include "option.h"
#include "input.h"
#include "rfpriv.h"
#include "hit.h"
#include "cluster.h"
#include "breakpoint.h"
#include "sample.h"
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <errno.h>
#include <zlib.h>
#include <set>



int main(int argc, char *argv[])
{
	if (argc < 2) {
		fprintf(stderr, "Usage: rf target.fa\n");
		return 1;
	}

	idxopt_t siopt(20, 1);
	idxopt_t liopt(100, 1);
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
	vector<sample> samples(genome.getSize());
	for (i = 0; i < genome.getSize(); ++i) {
		
		samples[i].process (i, i, mi_s, mi_l, fopts, siopt, liopt, 1);

		// // get diagonal cluster matches
		// vector<hit> dense_fhits, dense_rhits;
		// rf_hit(i, i, mi_s, dense_fhits, dense_rhits, fopts, siopt, 1); // find rhits also
		// cerr << "finish rf_hit for dense kmer!" << endl;
		// // clusters dense_clusts;
		// cleanDiag(dense_fhits, samples[i].dense_clusts, fopts, siopt, 0);
		// cleanDiag(dense_rhits, samples[i].dense_clusts, fopts, siopt, 1);
		// dense_fhits.clear();
		// dense_rhits.clear();
		// cerr << "finish cleanDiag and store diagonal clusters for sparse kmer!" << endl;

		// vector<hit> sparse_fhits, sparse_rhits;
		// rf_hit(i, i, mi_l, sparse_fhits, sparse_rhits, fopts, liopt, 1); // find rhits also
		// cerr << "finish rf_hit for dense kmer!" << endl;
		// // clusters sparse_clusts;
		// cleanDiag(sparse_fhits, samples[i].sparse_clusts, fopts, liopt, 0);
		// cleanDiag(sparse_rhits, samples[i].sparse_clusts, fopts, liopt, 1);
		// sparse_fhits.clear();
		// sparse_rhits.clear();
		// cerr << "finish cleanDiag and store diagonal clusters for sparse kmer!" << endl;

		// // flip the clusters
		// flipCluster(samples[i].dense_clusts, fopts, siopt, 1);
		// flipCluster(samples[i].sparse_clusts, fopts, liopt, 0);

		// // (TODO) Jingwen: make sure the code work for inversed cluster
		// // trim the breakpoints on y-axis
		// vector<uint32_t> bps_Y, bps_X;
		// trimOnY(samples[i].dense_clusts, samples[i].sparse_clusts, bps_Y, fopts);

		// // trim the breakpoints on x-axis;
		// trimOnX(samples[i].dense_clusts, samples[i].sparse_clusts, bps_Y, bps_X, fopts);
		// samples[i].breakpoints = bps_X;
		// bps_Y.clear();
		// bps_X.clear();






		// storeDiagCluster(fhits, clusts, 0, siopt, fopts);
		// storeDiagCluster(rhits, clusts, 1, siopt, fopts);
		// cerr << "finish storeDiagCluster!" << endl;

		// sort(clusts.begin(), clusts.end(), clustDiagonalSortOp);

		// // find overlap clusters && trim the edge of clusters && add back the main diagonal 
		// trim_ovpClusters(clusts, fopts.clusterTrimedge);
		// cerr << "finish trim_ovpClusters!" << endl;
		// // clusts.push_back(cluster(0, 0, 0, genoms->getLen(i), 0, genoms->getLen(i), 0));

		// // split clusters by endpoints
		// clusters splitclusts;
		// splitClusters(clusts, splitclusts);
		// trim_ovpClusters(splitclusts, 100);
		// cerr << "finish splitClusters!" << endl;

		// // get the fragment label for each sample
		// fragLabel(splitclusts, fopts);
		// cerr << "finish fragLabel!" << endl;






		// // find copy boundary for each sample based on the hits
		// clusters clusts;
		// cleanDiag(fhits);
		// storeDiagCluster(fhits, clusts, siopt, fopt);
		// sort(clusts.begin(), clusts.end(), clustDiagonalSortOp);
		// if (clusts.size == 0)
		// 	// do something 
		// int64_t mergeDiag = min((int64_t)clusts[0].xStart - (int64_t)clusts[0].xEnd - 1000, 3000);
		// clusters reclusts;
		// mergeDiagCluster(clusts, reclusts, mergeDiag, 0);

		// findCopyLoc(clusts); // copylen
		// findDel(clusts, reclusts);
		// fhits.clear();
		// rhits.clear();
	}


	// find large hits for each sample
	for (i = 0; i < genome.getSize(); ++i) {
		// check each small grid, unify them if a full cluster present; fragment them if a half cluster present; 
	}	


	// find large hits for a pair of sample
	for (i = 0; i < genome.getSize(); ++i) {
		// for each fragment in sample i, find the corresponding fragment in sample j; 
		// unify them if a full cluster present; fragment them if a half cluster present;
	}	
	return 0;
	










	// // open query file for reading; you may use your favorite FASTA parser
	// gzFile f = gzopen(argv[1], "r");
	// assert(f);
	// kseq_t *ks = kseq_init(f);



	// // generate two sets of index
	// mm_idx_reader_read(r_s, nthreads);
	// mm_idx_reader_read(r_d, nthreads);

	// // Read seq and collects hits
	// while (kseq_read(ks) >= 0) { // each kseq_read() call reads one query sequence
	// 	collect_hits();
	// 	repeat_frag();

	// }

	// while ((mi = mm_idx_reader_read(r, nthreads)) != 0) { // traverse each part of the index
	// 	// mm_mapopt_update(&mopt, mi); // this sets the maximum minimizer occurrence; TODO: set a better default in mm_mapopt_init()!
	// 	mm_tbuf_t *tbuf = mm_tbuf_init(); // thread buffer; for multi-threading, allocate one tbuf for each thread
	// 	gzrewind(f);
	// 	kseq_rewind(ks);
	// 	while (kseq_read(ks) >= 0) { // each kseq_read() call reads one query sequence
	// 		mm_reg1_t *reg;
	// 		int j, i, n_reg;
	// 		reg = mm_map(mi, ks->seq.l, ks->seq.s, &n_reg, tbuf, &mopt, 0); // get all hits for the query
	// 		for (j = 0; j < n_reg; ++j) { // traverse hits and print them out
	// 			mm_reg1_t *r = &reg[j];
	// 			assert(r->p); // with MM_F_CIGAR, this should not be NULL
	// 			printf("%s\t%d\t%d\t%d\t%c\t", ks->name.s, ks->seq.l, r->qs, r->qe, "+-"[r->rev]);
	// 			printf("%s\t%d\t%d\t%d\t%d\t%d\t%d\tcg:Z:", mi->seq[r->rid].name, mi->seq[r->rid].len, r->rs, r->re, r->mlen, r->blen, r->mapq);
	// 			for (i = 0; i < r->p->n_cigar; ++i) // IMPORTANT: this gives the CIGAR in the aligned regions. NO soft/hard clippings!
	// 				printf("%d%c", r->p->cigar[i]>>4, MM_CIGAR_STR[r->p->cigar[i]&0xf]);
	// 			putchar('\n');
	// 			free(r->p);
	// 		}
	// 		free(reg);
	// 	}
	// 	mm_tbuf_destroy(tbuf);
	// 	mm_idx_destroy(mi);
	// }
	// mm_idx_reader_close(r); // close the index reader
	// kseq_destroy(ks); // close the query file
	// gzclose(f);
	// return 0;
}
