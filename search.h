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

class searcher {
public:
	searcher();
	virtual ~searcher();
	void search8(BYTE * * q_start, BYTE gap_open_penalty, BYTE gap_extend_penalty, BYTE * score_matrix, BYTE * dprofile, BYTE * hearray,
			unsigned long sequences, unsigned long * seqnos, unsigned long * scores, unsigned long * diffs,
			unsigned long * alignmentlengths, unsigned long qlen, unsigned long dirbuffersize, unsigned long * dirbuffer, Db_data* db);

	void search16(WORD * * q_start, WORD gap_open_penalty, WORD gap_extend_penalty, WORD * score_matrix, WORD * dprofile, WORD * hearray,
			unsigned long sequences, unsigned long * seqnos, unsigned long * scores, unsigned long * diffs,
			unsigned long * alignmentlengths, unsigned long qlen, unsigned long dirbuffersize, unsigned long * dirbuffer, Db_data* db);
};

#endif /* SEARCH_H_ */
