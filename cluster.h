/*
 * Bigmatrix.h
 *
 *  Created on: Dec 24, 2014
 *      Author: mimitantono
 */

#ifndef CLUSTER_H_
#define CLUSTER_H_

#include <vector>
#include <map>
#include <pthread.h>
#include <queue>
#include "qgram.h"
#include "scan.h"
#include "db.h"
#include "clusterresult.h"
#include "util.h"
#include <locale.h>

class Cluster {
public:
	Cluster();
	virtual ~Cluster();
	void form_clusters();
	void print_debug();
	void print_clusters();
	void run_thread(int thread_id, int total_thread);
private:
	pthread_mutex_t row_id_mutex;
	pthread_mutex_t result_mutex;

	cluster_result result;
	class scanner * scanner;

	unsigned long int total_match;
	unsigned long int total_qgram;
	unsigned long int total_scan;
	unsigned long int row_full;
	unsigned long int row_reference;
	unsigned long int * row_stat;
	unsigned long int row_id_status;
	unsigned long int cluster_id;

	std::vector<unsigned long int> * targetampliconids;
	std::queue<unsigned long int> * next_step;
	std::queue<unsigned int> * next_step_level;
	std::vector<unsigned long int> ** next_comparison;
	std::map<unsigned long int, bool> * match_statistics;

	unsigned long int row_id_dispenser();
	void vector_put(int thread_id, unsigned long int row, unsigned long int col);
	void process_row(bool write_reference, bool use_reference, int thread_id, unsigned long int row_id, unsigned int iteration);
	void write_next_comparison(int thread_id, unsigned long int col, unsigned int distance);
};

#endif /* CLUSTER_H_ */
