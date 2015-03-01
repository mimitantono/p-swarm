#include "matrix.h"
#include <stdio.h>
#include <stdlib.h>
#include "util.h"
#include "property.h"

long SCORELIMIT_7 = 0;
long SCORELIMIT_8;
long SCORELIMIT_16;
long SCORELIMIT_32;
long SCORELIMIT_63;
char BIAS;

long * Matrix::score_matrix_63;
unsigned char * Matrix::score_matrix_8;
unsigned short * Matrix::score_matrix_16;
char sym_nt[] = "-acgt                           ";

void Matrix::score_matrix_dump() {
	fprintf(stderr, "     ");
	for (int i = 0; i < 16; i++)
		fprintf(stderr, "%2d", i);
	fprintf(stderr, "\n");
	fprintf(stderr, "     ");
	for (int i = 0; i < 16; i++)
		fprintf(stderr, " %c", sym_nt[i]);
	fprintf(stderr, "\n");
	for (int i = 0; i < 16; i++) {
		fprintf(stderr, "%2d %c ", i, sym_nt[i]);
		for (int j = 0; j < 16; j++) {
			fprintf(stderr, "%2ld", score_matrix_63[(i << 5) + j]);
		}
		fprintf(stderr, "\n");
	}
}

void Matrix::score_matrix_read() {
	int a, b;
	long sc, lo, hi;

	score_matrix_8 = new unsigned char[32 * 32];
	score_matrix_16 = new unsigned short[32 * 32];
	score_matrix_63 = new long[32 * 32];

	hi = -1000;
	lo = 1000;

	for (a = 0; a < 16; a++)
		for (b = 0; b < 16; b++) {
			sc = ((a == b) && (a > 0) && (b > 0)) ? 0 : Property::penalty_mismatch;
			// sc = (a==b) ? matchscore : mismatchscore;
			if (sc < lo)
				lo = sc;
			if (sc > hi)
				hi = sc;
			score_matrix_63[(a << 5) + b] = sc;
		}

	SCORELIMIT_8 = 256 - hi;
	SCORELIMIT_16 = 65536 - hi;

	for (a = 0; a < 32; a++)
		for (b = 0; b < 32; b++) {
			sc = score_matrix_63[(a << 5) + b];
			score_matrix_8[(a << 5) + b] = (unsigned char) sc;
			score_matrix_16[(a << 5) + b] = (unsigned short) sc;
		}
}

void Matrix::score_matrix_init() {
	score_matrix_read();
//	score_matrix_dump();
}

void Matrix::score_matrix_free() {
	if (Matrix::score_matrix_8)
		delete[] Matrix::score_matrix_8;
	if (Matrix::score_matrix_16)
		delete[] Matrix::score_matrix_16;
	if (Matrix::score_matrix_63)
		delete[] Matrix::score_matrix_63;
}

