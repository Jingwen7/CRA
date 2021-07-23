#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <errno.h>
#include <zlib.h>
#include "genome.h"
#include "index.h"
#include "option.h"
#include "input.h"
#include "rfpriv.h"


int main(int argc, char *argv[])
{
	if (argc < 2) {
		fprintf(stderr, "Usage: rf target.fa\n");
		return 1;
	}

	idxopt_t siopt(17, 1);
	idxopt_t liopt(100, 1);
	// fragopt_t fopt;

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
		mi_s.idx_gen();
		mi_s.idx_sort();
		mi_s.storeIndex(sidxFile);
	}
	if (mi_l.readIndex(lidxFile) == 1) {
		mi_l.idx_gen();
		mi_l.idx_sort();
		mi_l.storeIndex(sidxFile);
	}

	// find small&dense hits for each sample
	uint32_t i, j;
	vector<hit> fhits, rhits;
	for (i = 0; i < genome.getSize(); ++i) {
		hit(i, i, mi_s, fhits, rhits);

		// find copy boundary for each sample based on the hits



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
