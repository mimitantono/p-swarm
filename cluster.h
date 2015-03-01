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
#include "clusterresult.h"

class Cluster {
public:
	Cluster();
	virtual ~Cluster();
	void find_and_add_singletons();
	void print_debug();
	void print_clusters();
	void run_thread(int thread_id, int total_thread);
private:
	pthread_mutex_t row_id_mutex;
	pthread_mutex_t result_mutex;

	cluster_result result;
	class scanner * scanner;

	//status variables
	unsigned long int matches_found;
	unsigned long int qgram_performed;
	unsigned long int scan_performed;
	unsigned long int row_full;
	unsigned long int row_reference;
	unsigned long int current_row_id;
	unsigned long int current_cluster_id;
	unsigned long int * row_stat_by_thread;
	unsigned long int * row_stat_by_iteration;

	//temp variables (flags and etc)
	std::vector<unsigned long int> * targetampliconids;
	std::queue<unsigned long int> * next_step;
	std::queue<unsigned int> * next_step_level;
	std::vector<unsigned long int> ** next_comparison;
	std::map<unsigned long int, bool> * match_statistics;
//	bool * row_visited;

	unsigned long int get_next_row_id(int thread_id);
	void add_match_to_cluster(int thread_id, unsigned long int row, unsigned long int col);
	void process_row(bool write_reference, bool use_reference, int thread_id, unsigned long int row_id, unsigned int iteration);
	inline void write_next_comparison(int thread_id, unsigned long int col, unsigned int distance);
	void reset_flags(int thread_id);
	void qgram_diff_full_row(unsigned long int row_id, int thread_id, bool write_reference);
	void walkthrough_row_by_reference(unsigned int iteration, int thread_id, unsigned long int row_id);
	void prepare_alignment(unsigned long int col_id, unsigned long int row_id, int thread_id);
};

#endif /* CLUSTER_H_ */
