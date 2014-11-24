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

struct ampliconinfo_s {
	unsigned ampliconid;
	unsigned diffestimate; /* lower bound estimate of dist from initial seed */
	unsigned swarmid;
	unsigned generation;
	unsigned radius; /* actual diff from initial seed */
};

class cluster_job {
private:
	Db_data * db;
	class scanner scanner;
public:
	cluster_job(Db_data * db);
	virtual ~cluster_job();
	cluster_result * algo_run(int threadid);
};

#endif /* CLUSTER_H_ */
