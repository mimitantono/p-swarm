/*
 * Bigmatrix.h
 *
 *  Created on: Dec 24, 2014
 *      Author: mimitantono
 */

#ifndef BIGMATRIX_H_
#define BIGMATRIX_H_

#include <vector>
#include <pthread.h>
#include "qgram.h"
#include "scan.h"
#include "db.h"
#include "clusterresult.h"
#include "util.h"

class Bigmatrix {
public:
	Bigmatrix(Db_data * db);
	virtual ~Bigmatrix();
	void form_clusters();
	void print_matrix();
	void print_clusters();
	void calculate_partition(int thread_id, int total_thread);
	void init_partition(int thread_id, int total_thread);
private:
	search_data * search_data;
	Db_data * db;
	int ** matrix;
	cluster_result result;
	void crawl_row(bool ** row_guestbook, unsigned long int cluster_id, unsigned long int row_id, int generation);
	bool has_match(int row_id);
};

#endif /* BIGMATRIX_H_ */
