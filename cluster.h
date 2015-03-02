/*
 * Bigmatrix.h
 *
 *  Created on: Dec 24, 2014
 *      Author: mimitantono
 */

#ifndef CLUSTER_H_
#define CLUSTER_H_

#include <pthread.h>
#include "clusterresult.h"
#include "clusterdata.h"

class Cluster {
public:
	Cluster();
	virtual ~Cluster();
	void find_and_add_singletons();
	void print_debug(cluster_data ** cluster_data);
	void print_clusters();
	void run_thread(cluster_data * cluster_data, int total_thread);
private:
	pthread_mutex_t row_id_mutex;
	pthread_mutex_t result_mutex;

	cluster_result result;
	unsigned long int current_row_id;
	unsigned long int current_cluster_id;

	unsigned long int get_next_row_id();
	void add_match_to_cluster(cluster_data * cluster_data, unsigned long int row, unsigned long int col);
	void process_row(bool write_reference, bool use_reference, cluster_data * cluster_data, unsigned long int row_id,
			unsigned int iteration);
	void qgram_diff_full_row(unsigned long int row_id, cluster_data * cluster_data, bool write_reference);
	void walkthrough_row_by_reference(unsigned int iteration, cluster_data * cluster_data, unsigned long int row_id);
	void prepare_alignment(unsigned long int col_id, unsigned long int row_id, cluster_data * cluster_data);
};

#endif /* CLUSTER_H_ */
