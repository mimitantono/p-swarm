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
#include "qgram.h"
#include "scan.h"
#include "db.h"
#include "clusterresult.h"
#include "util.h"

struct eqstr {
	bool operator()(unsigned long int s1, unsigned long int s2) const {
		return s1 == s2;
	}
};

class Cluster {
public:
	Cluster();
	virtual ~Cluster();
	void form_clusters();
	void print_debug();
	void print_clusters();
	void run_thread(int thread_id, int total_thread);
private:
	pthread_mutex_t workmutex;
	class scanner * scanner;
	unsigned long int row_id_dispenser();
	unsigned long int ** targetampliconids;
	unsigned long int total_match;
	unsigned long int total_qgram;
	unsigned long int total_scan;
	unsigned long int row_full;
	unsigned long int row_reference;
	unsigned long int * row_stat;
	unsigned long int row_id_status;
	std::vector<unsigned long int> * matrix_x;
	std::vector<unsigned long int> * matrix_y;
	std::vector<unsigned long int> * next_step;
	std::vector<unsigned long int> * next_comparison;
	cluster_result result;
	void vector_put(int thread_id, unsigned long int row, unsigned long int col);
	void process_row(bool write_reference, bool use_reference, int thread_id, unsigned long int row_id);
};

#endif /* CLUSTER_H_ */
