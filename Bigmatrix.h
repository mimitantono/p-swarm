/*
 * Bigmatrix.h
 *
 *  Created on: Dec 24, 2014
 *      Author: mimitantono
 */

#ifndef BIGMATRIX_H_
#define BIGMATRIX_H_

#include <vector>
#include <array>
#include <unordered_map>
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
//	struct search_data * search_data;
	Db_data * db;
	class scanner * scanner;
	unsigned long int total_match;
	unsigned long int total_skip;
	std::vector<unsigned long int *> * matrix;
//	unsigned long int ** guestbook;
	cluster_result result;
	void crawl_row(bool ** row_guestbook, unsigned long int cluster_id, unsigned long int row_id, int generation);
	bool has_match(unsigned long int row_id);
	bool vector_contains(std::vector<unsigned long int *> * vector, unsigned long int row, unsigned long int col);
	void vector_put(std::vector<unsigned long int *> * vector, unsigned long int row, unsigned long int col);
};

#endif /* BIGMATRIX_H_ */
