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
#include "clusterresult.h"
#include "db.h"

struct partition_info {
	unsigned long start;
	unsigned long end;
	unsigned long threadid;
};

struct ampliconinfo_s {
	unsigned ampliconid;
	unsigned diffestimate; /* lower bound estimate of dist from initial seed */
	unsigned swarmid;
	unsigned generation;
	unsigned radius; /* actual diff from initial seed */
};

cluster_result * algo_run(partition_info partition);

#endif /* CLUSTER_H_ */
