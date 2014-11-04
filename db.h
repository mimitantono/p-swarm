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

using namespace std;

struct seqinfo_s
{
  char * header;
  char * seq;
  unsigned int headerlen;
  unsigned int headeridlen;
  unsigned int seqlen;
  unsigned int abundance;
  unsigned int clusterid;
  unsigned int hdrhash;
  int abundance_start;
  int abundance_end;
};

typedef struct seqinfo_s seqinfo_t;

typedef unsigned char qgramvector_t[QGRAMVECTORBYTES];
extern qgramvector_t * qgrams;
extern unsigned long sequences;

inline unsigned char * db_getqgramvector(unsigned long seqno) {
	return (unsigned char*) (qgrams + seqno);
}

void db_read(string filename);

unsigned long db_getsequencecount();
unsigned long db_getnucleotidecount();

unsigned long db_getlongestheader();
unsigned long db_getlongestsequence();

seqinfo_t * db_getseqinfo(unsigned long seqno);

char * db_getsequence(unsigned long seqno);
unsigned long db_getsequencelen(unsigned long seqno);

void db_getsequenceandlength(unsigned long seqno,
                             char ** address,
                             long * length);

char * db_getheader(unsigned long seqno);
unsigned long db_getheaderlen(unsigned long seqno);

unsigned long db_getabundance(unsigned long seqno);

void db_showsequence(unsigned long seqno);
void db_showall();
void db_print_info();
void db_free();

void db_putseq(long seqno);

void db_qgrams_init();
void db_qgrams_done();

void fprint_id(FILE * stream, unsigned long x);
void fprint_id_noabundance(FILE * stream, unsigned long x);

#endif /* DB_H_ */
