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

#define QGRAMLENGTH 5
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

struct seqinfo_s {
	char * header;
	char * seq;
	unsigned int headerlen;
	unsigned int headeridlen;
	unsigned int seqlen;
	unsigned int abundance;
	unsigned long int clusterid;
	unsigned int hdrhash;
	int abundance_start;
	int abundance_end;
};

typedef struct seqinfo_s seqinfo_t;
typedef unsigned char qgramvector_t[QGRAMVECTORBYTES];

class Db_data {
private:
	void showseq(char * seq);
	std::vector<seqinfo_t> seqindex;
	void qgrams_init();
	bool process_line(long line);
	static bool detect_duplicates(std::vector<Db_data*>& db);
public:
	qgramvector_t * qgrams;
	unsigned long sequences;
	unsigned long nucleotides;
	unsigned long longest;
	int threadid;

	static char* read_file(std::vector<Db_data*> &db, char * datap);

	Db_data();
	virtual ~Db_data();
	void print_debug();
	unsigned char * get_qgram_vector(unsigned long seq_no);
	seqinfo_t * get_seqinfo(unsigned long seqno);
	queryinfo_t get_sequence_and_length(unsigned long seqno);
	void show_sequence(unsigned long seqno);
	void show_all();
	void print_info();
	void put_seq(long seqno);
};

#endif /* DB_H_ */

