/*
 * parallel.cpp
 *
 *  Created on: Nov 2, 2014
 *      Author: mimitantono
 */

#include "parallel.h"

cluster_result * Parallel::results;

Parallel::Parallel() {
	results = new cluster_result[Property::partition];
}

Parallel::~Parallel() {
}

typedef struct thread_data {
	unsigned long thread_id;
	Db_data * db;
} thread_data;

void *run_cluster(void *threadargs) {
	thread_data *my_data = (thread_data*) threadargs;
	class cluster_job *cluster_job = new class cluster_job(my_data->db);
	cluster_job->algo_run(my_data->thread_id, &Parallel::results[my_data->thread_id]);
	delete cluster_job;
	pthread_exit(NULL);
}

void Parallel::run() {
	pthread_t threads[Property::partition];
	pthread_attr_t attr;
	int rc;
	long t;
	void *status;
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	thread_data *thread_data_array = new thread_data[Property::partition];
	std::vector<Db_data*> db;
	char * datap = (char *) xmalloc(MEMCHUNK);
	datap = Db_data::read_file(db, datap);

	for (long i = 0; i < Property::partition; i++) {
		thread_data_array[i].thread_id = (unsigned long) i;
		thread_data_array[i].db = db[i];
		rc = pthread_create(&threads[i], &attr, run_cluster, (void *) &thread_data_array[i]);
		if (rc) {
			fprintf(stderr, "Error: unable to create thread, %d", rc);
			exit(-1);
		}
	}
	pthread_attr_destroy(&attr);
	for (t = 0; t < Property::partition; t++) {
		rc = pthread_join(threads[t], &status);
		if (rc) {
			fprintf(stderr, "ERROR; return code from pthread_join() is %d\n", rc);
			exit(-1);
		}
		printf("\nMain: completed join with thread %ld produced [%lu] clusters\n", t, results[t].clusters.size());
	}
	delete[] thread_data_array;
	thread_data_array = NULL;
	if (Property::partition > 1) {
		merger merger(&results, Property::partition);
		merger.merge_groups();
		merger.final_merge();
		merger.merge_result.print(Property::outfile);
	} else {
		results[0].print(Property::outfile);
	}
	fprintf(stderr, "\nfinished.....");
	for (int i = 0; i < db.size(); i++) {
		delete db[i];
	}
	if (datap)
		free(datap);
}

