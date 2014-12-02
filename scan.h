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
	void set_db(Db_data * db);
	void search_do(unsigned long query_no, unsigned long listlength, unsigned long * targets, unsigned long * scores, unsigned long * diffs,
			unsigned long * alignlengths, long bits);
	void search_begin();
	queryinfo_t query;
private:
	class searcher searcher;
	struct search_data * sd;
	unsigned long master_next;
	unsigned long master_length;
	unsigned long remainingchunks;
	unsigned long dirbufferbytes;
	unsigned long * master_targets;
	unsigned long * master_scores;
	unsigned long * master_diffs;
	unsigned long * master_alignlengths;
	int master_bits;
	Db_data * db;
	void search_alloc();
	void search_delete();
	void search_init();
	void search_chunk(long bits);
	int search_getwork(unsigned long * countref, unsigned long * firstref);
	void master_dump();
	void search_worker_core();
};

#endif /* SCAN_H_ */
