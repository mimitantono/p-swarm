/*
 * main.h
 *
 *  Created on: Oct 24, 2014
 *      Author: mimitantono
 */

#ifndef MAIN_H_
#define MAIN_H_

#include "db.h"
#include "matrix.h"
#include "property.h"
#include "cluster.h"
#include <getopt.h>
#include "parallel.h"
#include <sys/time.h>
#include "stdlib.h"
#include "Bigmatrix.h"


typedef struct thread_data {
	unsigned long thread_id;
	class Bigmatrix *bigmatrix;
	pthread_mutex_t workmutex;
} thread_data;

void destroy();
void run();
void args_init(int argc, char** argv);
void args_usage();
void calculate_matrix(class Bigmatrix *bigmatrix);

#endif /* MAIN_H_ */
