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
#include "scan.h"
#include "db.h"

struct queryinfo {
	unsigned long qno;
	long len;
	char * seq;
};

typedef struct queryinfo queryinfo_t;
extern queryinfo_t query;

typedef unsigned int UINT32;
typedef unsigned short WORD;
typedef unsigned char BYTE;
typedef BYTE VECTOR[16];

extern unsigned long longestdbsequence;

void search8(BYTE * * q_start, BYTE gap_open_penalty, BYTE gap_extend_penalty, BYTE * score_matrix, BYTE * dprofile, BYTE * hearray,
		unsigned long sequences, unsigned long * seqnos, unsigned long * scores, unsigned long * diffs, unsigned long * alignmentlengths,
		unsigned long qlen, unsigned long dirbuffersize, unsigned long * dirbuffer, Db_data* db);

void search16(WORD * * q_start, WORD gap_open_penalty, WORD gap_extend_penalty, WORD * score_matrix, WORD * dprofile, WORD * hearray,
		unsigned long sequences, unsigned long * seqnos, unsigned long * scores, unsigned long * diffs, unsigned long * alignmentlengths,
		unsigned long qlen, unsigned long dirbuffersize, unsigned long * dirbuffer, Db_data* db);

#endif /* SEARCH_H_ */
