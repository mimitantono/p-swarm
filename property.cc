/*
 * property.cc
 *
 *  Created on: Oct 26, 2014
 *      Author: mimitantono
 */

#include "property.h"

char * Property::queryname;
char * Property::matrixname;
long Property::gapopen;
long Property::gapextend;
long Property::penalty_gapopen;
long Property::penalty_gapextend;
long Property::penalty_mismatch = 18;
long Property::gapopenextend;
long Property::matchscore;
long Property::mismatchscore;
long Property::threads;
long Property::penalty_factor;
char * Property::databasename;
long Property::resolution;
FILE * Property::outfile = stdout;

void Property::init() {
	matchscore = DEFAULT_MATCHSCORE;
	mismatchscore = DEFAULT_MISMATCHSCORE;
	gapopen = DEFAULT_GAPOPEN;
	gapextend = DEFAULT_GAPEXTEND;

	penalty_mismatch = 2 * matchscore - 2 * mismatchscore;
	penalty_gapopen = 2 * gapopen;
	penalty_gapextend = 2 * matchscore + gapextend;

	penalty_factor = gcd(gcd(penalty_mismatch, penalty_gapopen),
			penalty_gapextend);

	penalty_mismatch /= penalty_factor;
	penalty_gapopen /= penalty_factor;
	penalty_gapextend /= penalty_factor;
}

void Property::print() {
	fprintf(stderr, "queryname : %s\n", queryname);
	fprintf(stderr, "matrixname : %s\n", matrixname);
	fprintf(stderr, "databasename : %s\n", databasename);
	fprintf(stderr, "gapopen : %ld\n", gapopen);
	fprintf(stderr, "gapextend : %ld\n", gapextend);
	fprintf(stderr, "penalty_gapopen : %ld\n", penalty_gapopen);
	fprintf(stderr, "penalty_gapextend : %ld\n", penalty_gapextend);
	fprintf(stderr, "penalty_mismatch : %ld\n", penalty_mismatch);
	fprintf(stderr, "gapopenextend : %ld\n", gapopenextend);
	fprintf(stderr, "matchscore : %ld\n", matchscore);
	fprintf(stderr, "mismatchscore : %ld\n", mismatchscore);
	fprintf(stderr, "threads : %ld\n", threads);
	fprintf(stderr, "penalty_factor : %ld\n", penalty_factor);
	fprintf(stderr, "resolution : %ld\n", resolution);
}

