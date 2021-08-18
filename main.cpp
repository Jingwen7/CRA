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
	uint32_t n = genome.getSize();

	// Get the breakpoints for each sample
	vector<sample> samples(genome.getSize());
	for (i = 0; i < n; ++i) {
		cerr << "process sample: " << *genome.getName(i) << " " << i << endl;
		samples[i].init(i, genome.getName(i));
		samples[i].process (samples, i, i, mi_s, mi_l, fopts, siopt, liopt, 1);
		samples[i].dump(genome.getName(i), fopts);
	}

	// get the clusters for across-sample 
	vector<vector<acrosample>> acrosamples(n);
	for (i = 0; i < n; ++i) {
		acrosamples[i].resize(n - i);
		for (j = i + 1; j < n; ++j) {
			acrosamples[i][j - i - 1].init(i, j, &samples[i], &samples[j]);
			acrosamples[i][j - i - 1].across_process(samples, i, j, mi_s, mi_l, fopts, siopt, liopt, 0);
			acrosamples[i][j - i - 1].decideStrand();
			// if (fopts.debug)
			// 	acrosamples[i][j - i - 1].dump(genome.getName(i), genome.getName(j), fopts, genome.getLen(j));
		}
	}	

	// unify the breakpoints
	for (i = 0; i < n; ++i) {
		vector<uint32_t> bps_i, bps_j, bps_i_single;
		set<uint32_t> bpset_i, bpset_j;
	    uint32_t k, l;
	    for (k = 0; k < samples[i].breakpoints.size(); ++k) {
	    	bpset_i.insert(samples[i].breakpoints[k]);
	    }

		// project samples[j]'s breakpoint to sample[i]
		// bpset_i will contain all the breakpoints
		// cluster bpset_i 
		// project back to other samples
		bps_i.clear();
		for (j = 0; j < n; ++j) {
			bps_i_single.clear();
			if (i == j) continue;
			if (j > i) {
				bpset_j.clear();
				bps_j.clear();
			    for (k = 0; k < samples[j].breakpoints.size(); ++k) {
			    	bpset_j.insert(samples[j].breakpoints[k]);
			    }
				// step 1: trim the y-axis (sample[j])
				firstTrim(samples[j], bpset_j, acrosamples[i][j - i - 1].dense_clusts, acrosamples[i][j - i - 1].sparse_clusts, bps_j, fopts, 0, 1);
				samples[j].breakpoints = bps_j;

				if (fopts.debug) {
					checkBreakpoints_Clusters(bps_j, acrosamples[i][j - i - 1].dense_clusts, 0);
					checkBreakpoints_Clusters(bps_j, acrosamples[i][j - i - 1].sparse_clusts, 0);
				}

				// step 2: project breakpoints from sample[j] to sample[i]
				secondTrim (samples[i], acrosamples[i][j - i - 1].dense_clusts, acrosamples[i][j - i - 1].sparse_clusts, bps_j, bps_i_single, fopts, 0, 1);
			}
			else { // j < i
				bpset_j.clear();
				bps_j.clear();	
			    for (k = 0; k < samples[j].breakpoints.size(); ++k) {
			    	bpset_j.insert(samples[j].breakpoints[k]);
			    }			
			    // step 1: trim the x-axis (sample[j])
				firstTrim(samples[j], bpset_j, acrosamples[j][i - j - 1].dense_clusts, acrosamples[j][i - j - 1].sparse_clusts, bps_j, fopts, 1, 1);
				samples[j].breakpoints = bps_j;

				// step 2: project breakpoints from sample[j] to sample[i]
				secondTrim (samples[i], acrosamples[j][i - j - 1].dense_clusts, acrosamples[j][i - j - 1].sparse_clusts, bps_j, bps_i_single, fopts, 1, 1);
			}

			for (k = 0; k < bps_i_single.size(); ++k) {
				bps_i.push_back(bps_i_single[k]);
				bpset_i.insert(bps_i_single[k]);
			}

			if (j > i) {
				insertPoint (acrosamples[i][j - i - 1].dense_clusts, bpset_i, 0);
				insertPoint (acrosamples[i][j - i - 1].sparse_clusts, bpset_i, 0);
				insertPoint (acrosamples[i][j - i - 1].dense_clusts, bps_i, 0);
				insertPoint (acrosamples[i][j - i - 1].sparse_clusts, bps_i, 0);
			}
			else {
				insertPoint (acrosamples[j][i - j - 1].dense_clusts, bpset_i, 1);
				insertPoint (acrosamples[j][i - j - 1].sparse_clusts, bpset_i, 1);
				insertPoint (acrosamples[j][i - j - 1].dense_clusts, bps_i, 1);
				insertPoint (acrosamples[j][i - j - 1].sparse_clusts, bps_i, 1);	
			}
				
			if (fopts.debug) {
				if (j > i) { // samples[i] is x-axis
					checkBreakpoints_Clusters(bpset_i, acrosamples[i][j - i - 1].dense_clusts, 1);
					checkBreakpoints_Clusters(bpset_i, acrosamples[i][j - i - 1].sparse_clusts, 1);
				} 
				else {  // samples[j] is y-axis
					checkBreakpoints_Clusters(bpset_i, acrosamples[j][i - j - 1].dense_clusts, 0);
					checkBreakpoints_Clusters(bpset_i, acrosamples[j][i - j - 1].sparse_clusts, 0);					
				}
			}
		}

		// cluster breakpoints in bpset_i
		// trim dense_clusts, sparse_clusts in across-sample dotplot
		vector<uint32_t> trimInfo(bps_i.size());
		sort(bps_i.begin(), bps_i.end());
		iota(trimInfo.begin(), trimInfo.end(), 0);
		clusterBreakpoints(bps_i, trimInfo, fopts);
		vector<bool> remove(bps_i.size(), 0);

		for (j = 0; j < n; ++j) {
			if (i == j) continue;
			if (j > i) {
				trimClusters_nomodifybreakpoints(samples[i], bps_i, trimInfo, acrosamples[i][j - i - 1].dense_clusts, acrosamples[i][j - i - 1].sparse_clusts, remove, 1, 1);
			}
			else { 
				trimClusters_nomodifybreakpoints(samples[i], bps_i, trimInfo, acrosamples[j][i - j - 1].dense_clusts, acrosamples[j][i - j - 1].sparse_clusts, remove, 0, 1);
			}
		}
		REsize(bps_i, remove);
		samples[i].breakpoints = bps_i;

		if (fopts.debug) {
			for (j = 0; j < n; ++j) {
				if (i == j) continue;
				if (j > i) {
					checkBreakpoints_Clusters(bps_i, acrosamples[i][j - i - 1].dense_clusts, 1);
					checkBreakpoints_Clusters(bps_i, acrosamples[i][j - i - 1].sparse_clusts, 1);					
				}
				else {
					checkBreakpoints_Clusters(bps_i, acrosamples[j][i - j - 1].dense_clusts, 0);
					checkBreakpoints_Clusters(bps_i, acrosamples[j][i - j - 1].sparse_clusts, 0);					
				}
			}
		}

		if (fopts.debug) {
			ofstream fclust("trimlines_unify.bed", ios_base::app);
		  	for (auto& it : bps_i)
		    	fclust << it << "\t" << *(samples[i].readname) << endl;
			fclust.close();					
		}

		// project samples[i]'s breakpoints to all other samples
		vector<uint32_t> bps_j_proj;
		for (j = 0; j < n; ++j) {
			bps_j.clear();
			if (i == j) continue;
			bps_j_proj = samples[j].breakpoints;
			if (j > i) { // samples[j] is y-axis
				insertPoint (acrosamples[i][j - i - 1].dense_clusts, bps_j, 1);
				insertPoint (acrosamples[i][j - i - 1].sparse_clusts, bps_j, 1);
			}
			else { // samples[j] is x-axis
				insertPoint (acrosamples[j][i - j - 1].dense_clusts, bps_j, 0);
				insertPoint (acrosamples[j][i - j - 1].sparse_clusts, bps_j, 0);
			}			

			if (j > i) 
				secondTrim (samples[j], acrosamples[i][j - i - 1].dense_clusts, acrosamples[i][j - i - 1].sparse_clusts, bps_i, bps_j_proj, fopts, 1, 1);
			else 
				secondTrim (samples[j], acrosamples[j][i - j - 1].dense_clusts, acrosamples[j][i - j - 1].sparse_clusts, bps_i, bps_j_proj, fopts, 0, 1);


		    // cluster breakpoings around cluster's boundaries
			// bpset_j: cluster's boundaries
			// bpset_j_proj : projected breakpoints from samples[i] + samples[j].breakpoints: samples[j]'s original breakpoints
			updateBreakpointsBasedOnPivot(samples[j].breakpoints, bps_j, bps_j_proj, fopts);
			cerr << "finish!" << endl;
		}



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
