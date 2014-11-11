/*
 * main.cc
 *
 *  Created on: Oct 19, 2014
 *      Author: mimitantono
 */
#include "main.h"

using namespace std;

int main(int argc, char** argv) {
	Matrix::score_matrix_init();
	CPU_Info::cpu_features_detect();
	CPU_Info::cpu_features_show();
	Property::init();
	args_init(argc, argv);
	run();
	destroy();
}

void args_init(int argc, char **argv) {
	/* Set defaults */
	opterr = 1;

	char short_options[] = "d:ho:t:vm:p:g:e:s:u:braz";

	static struct option long_options[] = { { "differences", required_argument, NULL, 'd' }, { "help", no_argument, NULL, 'h' }, {
			"output-file", required_argument, NULL, 'o' }, { "threads", required_argument, NULL, 't' }, { "match-reward", required_argument,
	NULL, 'm' }, { "mismatch-penalty", required_argument, NULL, 'p' }, { "gap-opening-penalty", required_argument, NULL, 'g' }, {
			"gap-extension-penalty", required_argument, NULL, 'e' }, { 0, 0, 0, 0 } };

	int option_index = 0;
	int c;

	while ((c = getopt_long(argc, argv, short_options, long_options, &option_index)) != -1) {
		switch (c) {
		case 'd':
			/* differences (resolution) */
			Property::set_resolution(atol(optarg));
			break;
		case 'o':
			/* output-file */
			Property::set_outfile(optarg);
			break;
		case 't':
			/* threads */
			Property::set_threads(atol(optarg));
			break;
		case 'm':
			/* match reward */
			Property::set_matchscore(atol(optarg));
			break;
		case 'p':
			/* mismatch penalty */
			Property::set_mismatchscore(-atol(optarg));
			break;
			/* gap opening penalty */
		case 'g':
			Property::set_gapopen(atol(optarg));
			break;
		case 'e':
			/* gap extension penalty */
			Property::set_gapextend(atol(optarg));
			break;
		case 'h':
			/* help */
		default:
			args_usage();
			exit(1);
			break;
		}
	}

	if (optind < argc)
		Property::databasename = argv[optind];

	if (!Property::outfile)
		Property::set_outfile("result.log");
	Property::print();
}

void run() {
	Parallel parallel;
	parallel.run();
	for (int i = 0; i < Parallel::results.size(); i++) {
		Parallel::results[i]->print();
	}
	parallel.~Parallel();
}

void args_show() {
	fprintf(stderr, "Database file:     %s\n", Property::databasename.c_str());
	fprintf(stderr, "Output file:       %s\n", Property::outfilename.c_str());
	fprintf(stderr, "Resolution (d):    %ld\n", Property::resolution);
	fprintf(stderr, "Threads:           %ld\n", Property::threads);
	fprintf(stderr, "Scores:            match: %ld, mismatch: %ld\n", Property::matchscore, Property::mismatchscore);
	fprintf(stderr, "Gap penalties:     opening: %ld, extension: %ld\n", Property::gapopen, Property::gapextend);
	fprintf(stderr, "Converted costs:   mismatch: %ld, gap opening: %ld, gap extension: %ld\n", Property::penalty_mismatch,
			Property::penalty_gapopen, Property::penalty_gapextend);
}

void args_usage() {
	/*               0         1         2         3         4         5         6         7          */
	/*               01234567890123456789012345678901234567890123456789012345678901234567890123456789 */

	fprintf(stderr, "Usage: main [OPTIONS] [filename]\n");
	fprintf(stderr, "  -d, --differences INTEGER           resolution (1)\n");
	fprintf(stderr, "  -h, --help                          display this help and exit\n");
	fprintf(stderr, "  -o, --output-file FILENAME          output result filename (stdout)\n");
	fprintf(stderr, "  -t, --threads INTEGER               number of threads to use (1)\n");
	fprintf(stderr, "  -m, --match-reward INTEGER          reward for nucleotide match (5)\n");
	fprintf(stderr, "  -p, --mismatch-penalty INTEGER      penalty for nucleotide mismatch (4)\n");
	fprintf(stderr, "  -g, --gap-opening-penalty INTEGER   gap open penalty (12)\n");
	fprintf(stderr, "  -e, --gap-extension-penalty INTEGER gap extension penalty (4)\n");
	fprintf(stderr, "\n");
}

void destroy() {
	Matrix::score_matrix_free();
}
