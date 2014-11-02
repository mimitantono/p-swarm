/*
 * cluster.h
 *
 *  Created on: Oct 26, 2014
 *      Author: mimitantono
 */

#ifndef CLUSTER_H_
#define CLUSTER_H_

#include <stdio.h>
#include "qgram.h"
#include "scan.h"

typedef struct partition_info {
	unsigned long start;
	unsigned long end;
	unsigned long threadid;
} partition;

void algo_run(partition_info partition);

#endif /* CLUSTER_H_ */
