/*
 * scan.h
 *
 *  Created on: Oct 26, 2014
 *      Author: mimitantono
 */

#ifndef SCAN_H_
#define SCAN_H_

#include <pthread.h>
#include "property.h"
#include "util.h"
#include "matrix.h"
#include "search.h"
#include "db.h"

class scanner {
public:
	scanner();
	virtual ~scanner();
	void search_do(unsigned long query_no, std::vector<unsigned long int>* targets);
	void search_begin();
	queryinfo_t query;
	std::vector<search_result> master_result;
	struct search_data * sd;
private:
	class searcher searcher;
	unsigned long dirbufferbytes;
	void search_init();
	void search_chunk(std::vector<unsigned long int>* targets);
	int search_getwork(unsigned long * countref, unsigned long * firstref);
	void search_worker_core();
};

#endif /* SCAN_H_ */
