/*
 * search.h
 *
 *  Created on: Oct 26, 2014
 *      Author: mimitantono
 */

#ifndef SEARCH_H_
#define SEARCH_H_

#include "cpu_info.h"
#include "tmmintrin.h"
#include "util.h"
#include "db.h"

struct search_data {
	BYTE ** qtable;
	WORD ** qtable_w;

	BYTE * dprofile;
	WORD * dprofile_w;

	BYTE * hearray;

	unsigned long * dir_array;

	unsigned long target_count;
	unsigned long target_index;
};

class searcher {
public:
	searcher();
	virtual ~searcher();
	void search8(struct search_data *sd, unsigned long * seqnos, unsigned long * scores, unsigned long * diffs,
			unsigned long * alignmentlengths, queryinfo_t * query, unsigned long dirbuffersize, Db_data* db);

	void search16(struct search_data *sd, unsigned long * seqnos, unsigned long * scores, unsigned long * diffs,
			unsigned long * alignmentlengths, queryinfo_t * query, unsigned long dirbuffersize, Db_data* db);
};

#endif /* SEARCH_H_ */
