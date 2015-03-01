/*
 * main.h
 *
 *  Created on: Oct 24, 2014
 *      Author: mimitantono
 */

#ifndef MAIN_H_
#define MAIN_H_

#include "cluster.h"
#include "clusterdata.h"

typedef struct thread_data {
	unsigned long thread_id;
	class cluster_data * cluster_data;
	class Cluster *cluster;
} thread_data;

void destroy();
void run();
void args_init(int argc, char** argv);
void args_usage();
void calculate_matrix(class Cluster *cluster, class cluster_data ** cluster_data);

#endif /* MAIN_H_ */
