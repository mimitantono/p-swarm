/*
 * search.h
 *
 *  Created on: Oct 26, 2014
 *      Author: mimitantono
 */

#ifndef SEARCH_H_
#define SEARCH_H_

#include "util.h"
#include<vector>

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
	static void search8(struct search_data *sd, std::vector<unsigned long int> * targets, unsigned long **result,
			queryinfo_t * query, unsigned long dirbuffersize, long longest);

	static void search16(struct search_data *sd, std::vector<unsigned long int> * targets, unsigned long **result,
			queryinfo_t * query, unsigned long dirbuffersize, long longest);
};

#endif /* SEARCH_H_ */
