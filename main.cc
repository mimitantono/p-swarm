/*
 * main.cc
 *
 *  Created on: Oct 19, 2014
 *      Author: mimitantono
 */
#include "main.h"
#include "cpu_info.h"
#include "db.h"
#include "matrix.h"
#include "property.h"
#include <getopt.h>
#include <sys/time.h>
#include "stdlib.h"
#include "property.h"

using namespace std;

int main(int argc, char** argv) {
	dup2(1, 2);
	Property::init();
	Matrix::score_matrix_init();
	CPU_Info::cpu_features_detect();
	CPU_Info::cpu_features_show();
	args_init(argc, argv);
	run();
	destroy();
}

void args_init(int argc, char **argv) {
	/* Set defaults */
	opterr = 1;

	char short_options[] = "d:ho:t:y:vm:p:g:e:s:u:bra";

	static struct option long_options[] = {
											{ "differences", required_argument, NULL, 'd' },
											{ "help", no_argument, NULL, 'h' },
											{ "output-file", required_argument, NULL, 'o' },
											{ "threads", required_argument, NULL, 't' },
											{ "match-reward", required_argument, NULL, 'm' },
											{ "mismatch-penalty", required_argument, NULL, 'p' },
											{ "gap-opening-penalty", required_argument, NULL, 'g' },
											{ "gap-extension-penalty", required_argument, NULL, 'e' },
											{ "depth", required_argument, NULL, 'y' },
											{ 0, 0, 0, 0 }, };

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
		case 'y':
			Property::depth = atol(optarg);
			Property::enable_flag = true;
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
	Property::recalculate();
	Property::print();
}

void run() {
	timeval start, end;
	gettimeofday(&start, NULL);
	Property::db_data.read_file();
	class Cluster cluster;
	cluster_data * cluster_data = new class cluster_data[Property::threads];
	calculate_matrix(&cluster, &cluster_data);
	gettimeofday(&end, NULL);
	double dif1 = (end.tv_sec - start.tv_sec) * 1000 + (end.tv_usec - start.tv_usec) / 1000;
	gettimeofday(&start, NULL);
	cluster.print_debug(&cluster_data);
	cluster.find_and_add_singletons();
	cluster.print_clusters();
	gettimeofday(&end, NULL);
	double dif2 = (end.tv_sec - start.tv_sec) * 1000 + (end.tv_usec - start.tv_usec) / 1000;
	fprintf(stderr, "Time calculation   : %.2lf\n", dif1 / 1000);
	fprintf(stderr, "Time singletons    : %.2lf\n", dif2 / 1000);
	fprintf(stderr, "Total time         : %.2lf\n", (dif1 + dif2) / 1000);
	fprintf(Property::outfile, "\nCalculate matrix duration %.2lf secs", dif1);
	fprintf(Property::outfile, "\nForm cluster duration %.2lf secs\n", dif2);
	delete[] cluster_data;
}

void *run_thread(void *threadargs) {
	thread_data *my_data = (thread_data*) threadargs;
	my_data->cluster->run_thread(my_data->cluster_data, Property::threads);
	pthread_exit(NULL);
}

void calculate_matrix(class Cluster *cluster, class cluster_data ** cluster_data) {
	pthread_attr_t attr;
	int rc;
	void *status;
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	thread_data *thread_data_array = new thread_data[Property::threads];

	pthread_t threads[Property::threads];
	for (long i = 0; i < Property::threads; i++) {
		thread_data_array[i].thread_id = (unsigned long) i;
		thread_data_array[i].cluster = cluster;
		thread_data_array[i].cluster_data = &((*cluster_data)[i]);
		thread_data_array[i].cluster_data->thread_id = i;
		rc = pthread_create(&threads[i], &attr, run_thread, (void *) &thread_data_array[i]);
		if (rc) {
			fprintf(stderr, "Error: unable to create thread, %d", rc);
			exit(-1);
		}
	}
	for (long t = 0; t < Property::threads; t++) {
		rc = pthread_join(threads[t], &status);
		if (rc) {
			fprintf(stderr, "ERROR; return code from pthread_join() is %d\n", rc);
			exit(-1);
		}
	}
	pthread_attr_destroy(&attr);
	if (thread_data_array)
		delete[] thread_data_array;
	thread_data_array = NULL;
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
	fprintf(stderr, "  -o, --output-file FILENAME          output result filename (result.log)\n");
	fprintf(stderr, "  -t, --threads INTEGER               number of threads to use (1)\n");
	fprintf(stderr, "  -m, --match-reward INTEGER          reward for nucleotide match (5)\n");
	fprintf(stderr, "  -p, --mismatch-penalty INTEGER      penalty for nucleotide mismatch (4)\n");
	fprintf(stderr, "  -g, --gap-opening-penalty INTEGER   gap open penalty (12)\n");
	fprintf(stderr, "  -e, --gap-extension-penalty INTEGER gap extension penalty (4)\n");
	fprintf(stderr, "  -y, --depth INTEGER                 max level for economic search\n");
	fprintf(stderr, "\n");
}

void destroy() {
	Matrix::score_matrix_free();
}
