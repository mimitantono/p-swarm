/*
 * property.cc
 *
 *  Created on: Oct 26, 2014
 *      Author: mimitantono
 */

#include "property.h"

long Property::gapopen;
long Property::gapextend;
long Property::penalty_gapopen;
long Property::penalty_gapextend;
long Property::penalty_mismatch = 18;
long Property::matchscore;
long Property::mismatchscore;
long Property::partition;
long Property::threads;
long Property::penalty_factor;
unsigned long Property::longest;
unsigned long Property::resolution;
unsigned long Property::diff_saturation;
unsigned long Property::bits;
FILE * Property::outfile;
FILE * Property::dbdebug;
std::string Property::databasename;
std::string Property::outfilename;

void Property::init() {
	matchscore = DEFAULT_MATCHSCORE;
	mismatchscore = DEFAULT_MISMATCHSCORE;
	gapopen = DEFAULT_GAPOPEN;
	gapextend = DEFAULT_GAPEXTEND;
	resolution = DEFAULT_RESOLUTION;
	partition = 1;
	threads = DEFAULT_THREADS;
	calculate_penalty();
	dbdebug = fopen("db.log", "w");
}

void Property::calculate_penalty() {
	if ((gapopen + gapextend) < 1)
		fatal("Illegal gap penalties specified");
	penalty_mismatch = 2 * matchscore - 2 * mismatchscore;
	penalty_gapopen = 2 * gapopen;
	penalty_gapextend = 2 * matchscore + gapextend;

	penalty_factor = gcd(gcd(penalty_mismatch, penalty_gapopen), penalty_gapextend);

	penalty_gapextend /= penalty_factor;
	penalty_mismatch /= penalty_factor;
	penalty_gapopen /= penalty_factor;
	diff_saturation = MIN(255 / penalty_mismatch, 255 / (penalty_gapopen + penalty_gapextend));
	if (Property::resolution <= diff_saturation)
		bits = 8;
	else
		bits = 16;
}

void Property::print() {
	fprintf(stderr, "databasename       : %s\n", databasename.c_str());
	fprintf(stderr, "outfilename        : %s\n", outfilename.c_str());
	fprintf(stderr, "gapopen            : %ld\n", gapopen);
	fprintf(stderr, "gapextend          : %ld\n", gapextend);
	fprintf(stderr, "matchscore         : %ld\n", matchscore);
	fprintf(stderr, "mismatchscore      : %ld\n", mismatchscore);
	fprintf(stderr, "penalty_gapopen    : %ld\n", penalty_gapopen);
	fprintf(stderr, "penalty_gapextend  : %ld\n", penalty_gapextend);
	fprintf(stderr, "penalty_mismatch   : %ld\n", Property::penalty_mismatch);
	fprintf(stderr, "penalty_factor     : %ld\n", penalty_factor);
	fprintf(stderr, "threads            : %ld\n", partition);
	fprintf(stderr, "resolution         : %ld\n", resolution);
}

void Property::set_resolution(long value) {
	resolution = value;
	if (resolution < 1)
		fatal("Illegal resolution specified");
}

void Property::set_outfile(std::string value) {
	outfilename = value;
	outfile = fopen(outfilename.c_str(), "w");
	if (!outfile)
		fatal("Unable to open output file for writing");

}
void Property::set_threads(long value) {
	threads = value;
	if ((threads < 1) || (threads > MAX_THREADS))
		fatal("Illegal number of threads specified");
}
void Property::set_matchscore(long value) {
	matchscore = value;
	if (matchscore < 1)
		fatal("Illegal match reward specified");
}
void Property::set_mismatchscore(long value) {
	mismatchscore = value;
	if (mismatchscore > -1)
		fatal("Illegal mismatch penalty specified");
}
void Property::set_gapopen(long value) {
	if (gapopen < 0)
		fatal("Illegal gap open specified.");
}
void Property::set_gapextend(long value) {
	if (gapextend < 0)
		fatal("Illegal gap extend specified.");
}

