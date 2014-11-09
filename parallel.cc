/*
 * parallel.cpp
 *
 *  Created on: Nov 2, 2014
 *      Author: mimitantono
 */

#include "parallel.h"

vector<cluster_result*> Parallel::results;

Parallel::Parallel() {
}

Parallel::~Parallel() {
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
	Parallel::results.push_back(algo_run(partition));
	pthread_exit(NULL);
}

void Parallel::run() {
	pthread_t threads[Property::threads];
	pthread_attr_t attr;
	int rc;
	long t;
	void *status;
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

	for (long i = 0; i < Property::threads; i++) {
		rc = pthread_create(&threads[i], &attr, run_cluster, (void *) (unsigned long) i);
		if (rc) {
			fprintf(stderr, "Error: unable to create thread, %d", rc);
			exit(-1);
		}
	}
	pthread_attr_destroy(&attr);
	for (t = 0; t < Property::threads; t++) {
		rc = pthread_join(threads[t], &status);
		if (rc) {
			fprintf(stderr, "ERROR; return code from pthread_join() is %d\n", rc);
			exit(-1);
		}
		printf("\nMain: completed join with thread %ld having a status of %ld\n", t + 1, (long) status);
	}
}

