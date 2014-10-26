/*
 * property.h
 *
 *  Created on: Oct 26, 2014
 *      Author: mimitantono
 */

#ifndef PROPERTY_H_
#define PROPERTY_H_

#include <stdio.h>
#include "util.h"

#define DEFAULT_GAPOPEN 12
#define DEFAULT_GAPEXTEND 4
#define DEFAULT_MATCHSCORE 5
#define DEFAULT_MISMATCHSCORE (-4)
#define DEFAULT_THREADS 1
#define DEFAULT_RESOLUTION 1
#define DEFAULT_ALTERNATIVE_ALGORITHM 0

class Property {
public:
	static void init();
	static void print();
	static char * queryname;
	static char * matrixname;
	static long gapopen;
	static long gapextend;
	static long penalty_gapopen;
	static long penalty_gapextend;
	static long penalty_mismatch;
	static long penalty_factor;
	static long gapopenextend;
	static long matchscore;
	static long mismatchscore;
	static long threads;
	static char * databasename;
	static long resolution;
	static FILE* outfile;
};



#endif /* PROPERTY_H_ */
