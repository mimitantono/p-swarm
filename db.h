/*
 * db.h
 *
 *  Created on: Oct 25, 2014
 *      Author: mimitantono
 */

#ifndef DB_H_
#define DB_H_

#define LINE_MAX 2048
#define HASH hash_cityhash64

#define MEMCHUNK 1048576
#define LINEALLOC LINE_MAX

#define QGRAMLENGTH 4
#define QGRAMVECTORBITS (1<<(2*QGRAMLENGTH))
#define QGRAMVECTORBYTES (QGRAMVECTORBITS/8)

#include "qgram.h"
#include "matrix.h"
#include "util.h"
#include <stdio.h>
#include <limits.h>
#include <regex.h>
#include <string>
#include <vector>

typedef struct seqinfo_s seqinfo_t;
typedef unsigned char qgramvector_t[QGRAMVECTORBYTES];

struct seqinfo_s {
	char * header;
	char * seq;
	unsigned int headerlen;
	unsigned int seqlen;
	unsigned int abundance;
	unsigned int hdrhash;
	qgramvector_t qgram;
};


class Db_data {
private:
	void showseq(char * seq);
	std::vector<seqinfo_t> seqindex;
	void qgrams_init();
	bool detect_duplicates();
	char * datap;
public:
	unsigned long sequences;
	unsigned long nucleotides;
	unsigned long longest;
	int threadid;

	void read_file();

	Db_data();
	virtual ~Db_data();
	void print_debug();
	seqinfo_t * get_seqinfo(unsigned long seqno);
	queryinfo_t get_queryinfo(unsigned long seqno);
	void show_sequence(unsigned long seqno);
	void show_all();
	void print_info();
	void put_seq(long seqno);
};

#endif /* DB_H_ */

