/*
 * parallel.cpp
 *
 *  Created on: Nov 2, 2014
 *      Author: mimitantono
 */

#include "parallel.h"

vector<cluster_result*> Parallel::results;

Parallel::Parallel(class Db_data * db) {
	this->db = db;
}

Parallel::~Parallel() {
}

long get_start(unsigned long threadid, Db_data * db) {
	return db->sequences / Property::threads * threadid;
}

long get_end(unsigned long threadid, Db_data* db) {
	if (threadid < Property::threads - 1) {
		return db->sequences / Property::threads * (threadid + 1) - 1;
	} else {
		return db->sequences - 1;
	}
}

typedef struct thread_data {
	unsigned long thread_id;
	Db_data * db;
} thread_data;

void *run_cluster(void *threadargs) {
	thread_data *my_data = (thread_data*) threadargs;
	fprintf(stderr, "\nMy Data: %lu", (unsigned long) my_data->thread_id);
	partition_info partition;
	partition.threadid = my_data->thread_id;
	partition.start = get_start(partition.threadid, my_data->db);
	partition.end = get_end(partition.threadid, my_data->db);
	fprintf(stderr, "\nStarting thread #%lu from %lu to %lu", partition.threadid, partition.start, partition.end);
	cluster_job cluster_job(my_data->db);
	Parallel::results.push_back(cluster_job.algo_run(partition));
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
	thread_data *thread_data_array = new thread_data[Property::threads];

	for (long i = 0; i < Property::threads; i++) {
		thread_data_array[i].thread_id = (unsigned long) i;
		thread_data_array[i].db = db;
		rc = pthread_create(&threads[i], &attr, run_cluster, (void *) &thread_data_array[i]);
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

