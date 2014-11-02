/*
 * parallel.cpp
 *
 *  Created on: Nov 2, 2014
 *      Author: mimitantono
 */

#include "parallel.h"

Parallel::Parallel() {
}

Parallel::~Parallel() {
	// TODO Auto-generated destructor stub
}

long get_start(unsigned long threadid) {
	return db_getsequencecount() / Property::threads * threadid;
}

long get_end(unsigned long threadid) {
	if (threadid < Property::threads - 1) {
		return db_getsequencecount() / Property::threads * (threadid + 1);
	} else {
		return db_getsequencecount();
	}
}

void *run_cluster(void *threadid) {
	partition_info partition;
	partition.threadid = (long) threadid;
	partition.start = get_start((long) threadid);
	partition.end = get_end((long) threadid);
	algo_run(partition);
	pthread_exit(NULL);
}

void Parallel::run() {
	pthread_t threads[Property::threads];
	int rc;
	for (long i = 0; i < Property::threads; i++) {
		rc = pthread_create(&threads[i], NULL, run_cluster, (void *) (unsigned long) i);
		if (rc) {
			fprintf(stderr, "Error: unable to create thread, %d", rc);
			exit(-1);
		}
	}
	pthread_exit(NULL);
}

