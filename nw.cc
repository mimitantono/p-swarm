#include "nw.h"

namespace AlignAttributes {
const unsigned char maskup = 1;
const unsigned char maskleft = 2;
const unsigned char maskextup = 4;
const unsigned char maskextleft = 8;
}

NW_Result::NW_Result() {
	nwscore = 0;
	nwdiff = 0;
	nwalignmentlength = 0;
	nwalignment = NULL;
}

NeedlemanWunsch::NeedlemanWunsch(unsigned long gapopen, unsigned long gapextend,
		int longestAmplicon) {
	this->gapopen = gapopen;
	this->gapextend = gapextend;
	this->dir = (unsigned char *) xmalloc(longestAmplicon * longestAmplicon);
	this->hearray = (unsigned long *) xmalloc(
			2 * longestAmplicon * sizeof(unsigned long));
}

//	Needleman/Wunsch/Sellers aligner
//
//		 finds a global alignment with minimum cost
//		 there should be positive costs/penalties for gaps and for mismatches
//		 matches should have zero cost (0)
//
//		 alignment priority when backtracking (from lower right corner):
//		 1. left/insert/e (gap in query sequence (qseq))
//		 2. align/diag/h (match/mismatch)
//		 3. up/delete/f (gap in database sequence (dseq))
//
//		 qseq: the reference/query/upper/vertical/from sequence
//		 dseq: the sample/database/lower/horisontal/to sequence
//
//		 typical costs:
//		 match: 0
//		 mismatch: 3
//		 gapopen: 4
//		 gapextend: 3
//
//		 input
//
//		 dseq: pointer to start of database sequence
//		 dend: pointer after database sequence
//		 qseq: pointer to start of query sequence
//		 qend: pointer after database sequence
//		 score_matrix: 32x32 matrix of longs with scores for aligning two symbols
//		 gapopen: positive number indicating penalty for opening a gap of length zero
//		 gapextend: positive number indicating penalty for extending a gap
//
//		 output
//
//		 nwscore: the global alignment score
//		 nwdiff: number of non-identical nucleotides in one optimal global alignment
//		 nwalignmentlength: the length of one optimal alignment
//		 nwalignment: cigar string with one optimal alignment
void NeedlemanWunsch::align(char * dseq, char * qseq, NW_Result * nw_result) {
	/* dir must point to at least qlen*dlen bytes of allocated memory
	 hearray must point to at least 2*qlen longs of allocated memory (8*qlen bytes) */

	long h, n, e, f;
	long unsigned *hep;

	long qlen = strlen(qseq);
	long dlen = strlen(dseq);

	memset(dir, 0, qlen * dlen);

	long i, j;

	for (i = 0; i < qlen; i++) {
		hearray[2 * i] = 1 * gapopen + (i + 1) * gapextend; // H (N)
		hearray[2 * i + 1] = 2 * gapopen + (i + 2) * gapextend; // E
	}

	for (j = 0; j < dlen; j++) {
		hep = hearray;
		f = 2 * gapopen + (j + 2) * gapextend;
		h = (j == 0) ? 0 : (gapopen + j * gapextend);

		for (i = 0; i < qlen; i++) {
			long index = qlen * j + i;

			n = *hep;
			e = *(hep + 1);
			h += Matrix::score_matrix_63[(dseq[j] << 5) + qseq[i]];

			dir[index] |= (f < h ? AlignAttributes::maskup : 0);
			h = MIN(h, f);
			h = MIN(h, e);
			dir[index] |= (e == h ? AlignAttributes::maskleft : 0);

			*hep = h;

			h += gapopen + gapextend;
			e += gapextend;
			f += gapextend;

			dir[index] |= (f < h ? AlignAttributes::maskextup : 0);
			dir[index] |= (e < h ? AlignAttributes::maskextleft : 0);
			f = MIN(h, f);
			e = MIN(h, e);

			*(hep + 1) = e;
			h = n;
			hep += 2;
		}
	}

	long dist = hearray[2 * qlen - 2];

	/* backtrack: count differences and save alignment in cigar string */

	long score = 0;
	long alength = 0;
	long matches = 0;

	char * cigar = (char *) xmalloc(qlen + dlen + 1);
	char * cigarend = cigar + qlen + dlen + 1;

	char op = 0;
	int count = 0;
	*(--cigarend) = 0;

	i = qlen;
	j = dlen;

	while ((i > 0) && (j > 0)) {
		int d = dir[qlen * (j - 1) + (i - 1)];

		alength++;

		if ((op == 'I') && (d & AlignAttributes::maskextleft)) {
			score += gapextend;
			j--;
			pushop('I', &cigarend, &op, &count);
		} else if ((op == 'D') && (d & AlignAttributes::maskextup)) {
			score += gapextend;
			i--;
			pushop('D', &cigarend, &op, &count);
		} else if (d & AlignAttributes::maskleft) {
			score += gapextend;
			if (op != 'I')
				score += gapopen;
			j--;
			pushop('I', &cigarend, &op, &count);
		} else if (d & AlignAttributes::maskup) {
			score += gapextend;
			if (op != 'D')
				score += gapopen;
			i--;
			pushop('D', &cigarend, &op, &count);
		} else {
			score += Matrix::score_matrix_63[(dseq[j - 1] << 5) + qseq[i - 1]];
			if (qseq[i - 1] == dseq[j - 1])
				matches++;
			i--;
			j--;
			pushop('M', &cigarend, &op, &count);
		}
	}

	while (i > 0) {
		alength++;
		score += gapextend;
		if (op != 'D')
			score += gapopen;
		i--;
		pushop('D', &cigarend, &op, &count);
	}

	while (j > 0) {
		alength++;
		score += gapextend;
		if (op != 'I')
			score += gapopen;
		j--;
		pushop('I', &cigarend, &op, &count);
	}

	finishop(&cigarend, &op, &count);

	/* move and reallocate cigar */

	long cigarlength = cigar + qlen + dlen - cigarend;
	memmove(cigar, cigarend, cigarlength + 1);
	cigar = (char*) realloc(cigar, cigarlength + 1);

	nw_result->nwscore = dist;
	nw_result->nwdiff = alength - matches;
	nw_result->nwalignmentlength = alength;
	nw_result->nwalignment = cigar;

	if (score != dist) {
//			fprintf(stderr, "Error with query no %lu and db sequence no %lu:\n",
//					queryno, dbseqno);
		fprintf(stderr,
				"Initial and recomputed alignment score disagreement: %lu %lu\n",
				dist, score);
		fprintf(stderr, "Alignment: %s\n", cigar);
	}
}

void NeedlemanWunsch::pushop(char newop, char ** cigarendp, char * op,
		int * count) {
	if (newop == *op)
		(*count)++;
	else {
		*--*cigarendp = *op;
		if (*count > 1) {
			char buf[25];
			int len = sprintf(buf, "%d", *count);
			*cigarendp -= len;
			strncpy(*cigarendp, buf, len);
		}
		*op = newop;
		*count = 1;
	}
}

void NeedlemanWunsch::finishop(char ** cigarendp, char * op, int * count) {
	if ((op) && (count)) {
		*--*cigarendp = *op;
		if (*count > 1) {
			char buf[25];
			int len = sprintf(buf, "%d", *count);
			*cigarendp -= len;
			strncpy(*cigarendp, buf, len);
		}
		*op = 0;
		*count = 0;
	}
}
