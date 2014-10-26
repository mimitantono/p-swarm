/*
 * nw.h
 *
 *  Created on: Oct 21, 2014
 *      Author: mimitantono
 */
#ifndef NW_H_
#define NW_H_

#include <stdio.h>
#include <string.h>
#include <pthread.h>
#include <getopt.h>
#include <tmmintrin.h>
#include <stdlib.h>
#include <regex.h>
#include <limits.h>
#include <stdint.h>
#include "matrix.h"
#include "util.h"

class NW_Result {
public:
	NW_Result();
	unsigned long nwdiff;
	unsigned long nwscore;
	unsigned long nwalignmentlength;
	char * nwalignment;
};

class NeedlemanWunsch {
private:
	unsigned long gapopen;
	unsigned long gapextend;
	unsigned char * dir;
	unsigned long * hearray;
	void pushop(char newop, char ** cigarendp, char * op, int * count);
	void finishop(char ** cigarendp, char * op, int * count);
public:
	void align(char * dseq, char * qseq, NW_Result * nw_result);
	NeedlemanWunsch(unsigned long gapopen, unsigned long gapextend,
			int longestAmplicon);
};

#endif /* NW_H_ */
