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

struct search_result {
	unsigned long score;
	unsigned long diff;
	unsigned long alignlength;
};

class searcher {
public:
	searcher();
	virtual ~searcher();
	static void search8(struct search_data *sd, std::vector<queryinfo_t> targets, std::vector<search_result> *result, queryinfo_t * query,
			unsigned long dirbuffersize, long longest);

	static void search16(struct search_data *sd, std::vector<queryinfo_t> targets, std::vector<search_result> *result, queryinfo_t * query,
			unsigned long dirbuffersize, long longest);
};

#endif /* SEARCH_H_ */
