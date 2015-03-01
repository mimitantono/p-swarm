/*
 * util.h
 *
 *  Created on: Oct 25, 2014
 *      Author: mimitantono
 */

#ifndef UTIL_H_
#define UTIL_H_

#include <stdlib.h>

#define MIN(x,y) ((x)<(y)?(x):(y))

typedef unsigned int UINT32;
typedef unsigned short WORD;
typedef unsigned char BYTE;
typedef BYTE VECTOR[16];

struct queryinfo_t {
	unsigned long qno;
	long len;
	char * seq;
};

long gcd(long a, long b);
void fatal(const char * msg);
void fatal(const char * format, const char * message);
void * xmalloc(size_t size);
void * xrealloc(void * ptr, size_t size);
char * xstrchrnul(char *s, int c);
unsigned long hash_fnv_1a_64(unsigned char * s, unsigned long n);
unsigned int hash_fnv_1a_32(unsigned char * s, unsigned long n);
unsigned long hash_djb2(unsigned char * s, unsigned long n);
unsigned long hash_djb2a(unsigned char * s, unsigned long n);
unsigned long hash_cityhash64(unsigned char * s, unsigned long n);
void progress_init(const char * prompt, unsigned long size);
void progress_update(unsigned long progress);
void progress_done();

#endif /* UTIL_H_ */
