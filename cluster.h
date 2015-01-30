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

class Bigmatrix {
public:
	Bigmatrix();
	virtual ~Bigmatrix();
	void form_clusters();
	void print_debug();
	void print_clusters();
	void calculate_partition(int thread_id, int total_thread);
private:
	pthread_mutex_t workmutex;
	class scanner * scanner;
	unsigned long int ** targetampliconids;
	unsigned long int total_match;
	unsigned long int total_qgram;
	unsigned long int total_scan;
	unsigned long int temp_written;
	unsigned long int temp_cleaned;
	unsigned long int total_data;
	std::vector<unsigned long int> * matrix_x;
	std::vector<unsigned long int> * matrix_y;
	std::map<unsigned long int, std::vector<unsigned long int> > next_comparison;
	std::map<unsigned long int, int> comparison_log;
	cluster_result result;
	void vector_put(int thread_id, unsigned long int row, unsigned long int col);
};

#endif /* CLUSTER_H_ */
