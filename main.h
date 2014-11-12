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

void destroy();
void run(Db_data * db);
void args_init(int argc, char** argv);
void args_usage();

#endif /* MAIN_H_ */
