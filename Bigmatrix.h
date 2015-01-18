/*
 * Bigmatrix.h
 *
 *  Created on: Dec 24, 2014
 *      Author: mimitantono
 */

#ifndef BIGMATRIX_H_
#define BIGMATRIX_H_

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
	Bigmatrix(Db_data * db);
	virtual ~Bigmatrix();
	void form_clusters();
	void print_debug();
	void print_clusters();
	void calculate_partition(int thread_id, int total_thread, pthread_mutex_t workmutex);
	void init_partition(int thread_id, int total_thread);
private:
//	struct search_data * search_data;
	Db_data * db;
	class scanner * scanner;
	unsigned long int total_match;
	unsigned long int total_skip;
	unsigned long int total_qgram;
	unsigned long int total_scan;
	unsigned long int skip_by_guestbook;
	std::vector<unsigned long int *> * matrix;
	unsigned long int * guestbook;
	cluster_result result;
	void vector_put(std::vector<unsigned long int *> * vector, unsigned long int row, unsigned long int col);
	unsigned long int get_combined(unsigned long int a , unsigned long int b);
	bool skip_rowcol(unsigned long int a, unsigned long int b);
};

#endif /* BIGMATRIX_H_ */
